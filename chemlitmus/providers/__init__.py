"""Multi-database compound lookup.

::

    from chemlitmus.providers import resolve, get_provider, PROVIDERS

    res = resolve("aspirin")                      # all providers in parallel
    res = resolve("CHEMBL25", sources=["chembl", "chebi"])
    rec = get_provider("kegg").lookup("C01405")

A :class:`ResolveResult` carries one :class:`CompoundRecord` per source that answered, a merged
view that prefers the most complete fields, cross-references from UniChem, and a note on
whether the sources agree about what compound the query refers to.
"""

from __future__ import annotations

import concurrent.futures as cf
import time
from collections import Counter
from typing import Dict, Iterable, List, Optional, Type

from pydantic import BaseModel, Field

from chemlitmus.logging_config import logger
from chemlitmus.providers.base import QUERY_TYPES, CompoundRecord, Provider, ProviderError, _looks_like_inchikey
from chemlitmus.providers.chebi import ChEBIProvider
from chemlitmus.providers.chembl import ChEMBLProvider
from chemlitmus.providers.kegg import KEGGProvider
from chemlitmus.providers.pubchem import PubChemProvider
from chemlitmus.providers.unichem import unichem_xrefs

PROVIDERS: Dict[str, Type[Provider]] = {
    PubChemProvider.key: PubChemProvider,
    ChEMBLProvider.key: ChEMBLProvider,
    ChEBIProvider.key: ChEBIProvider,
    KEGGProvider.key: KEGGProvider,
}
"""Registered providers, keyed by short name."""

DEFAULT_SOURCES: List[str] = ["pubchem", "chembl", "chebi"]
"""Sources queried when none are specified. KEGG is opt-in: it has no structure search and is
slower, so it is most useful when you already know you want pathway context."""

_instances: Dict[str, Provider] = {}


def get_provider(key: str) -> Provider:
    k = key.lower()
    if k not in PROVIDERS:
        raise ValueError(f"Unknown source {key!r}. Available: {sorted(PROVIDERS)}")
    if k not in _instances:
        _instances[k] = PROVIDERS[k]()
    return _instances[k]


class SourceOutcome(BaseModel):
    source: str
    status: str = Field(description="'found' | 'not found' | 'error'")
    error: Optional[str] = None
    seconds: float = 0.0


class ResolveResult(BaseModel):
    """Everything learned about one query across the requested sources."""

    query: str
    query_type: str
    sources_requested: List[str]
    outcomes: List[SourceOutcome]
    records: List[CompoundRecord] = Field(default_factory=list, description="One per source that answered.")
    merged: Optional[CompoundRecord] = Field(None, description="Field-wise merge preferring the first non-empty value in source order.")
    cross_refs: Dict[str, str] = Field(default_factory=dict, description="Identifiers across databases (from the records plus UniChem).")
    consensus_inchikey: Optional[str] = None
    agreement: str = Field("unknown", description="'agree' when every record shares the consensus InChIKey, 'disagree' otherwise, 'unknown' when fewer than two records carry one.")

    @property
    def found(self) -> bool:
        return bool(self.records)

    def by_source(self, key: str) -> Optional[CompoundRecord]:
        return next((r for r in self.records if r.source == key), None)


_MERGE_FIELDS = ("name", "smiles", "inchi", "inchikey", "formula", "molecular_weight", "monoisotopic_mass",
                 "charge", "xlogp", "hbd", "hba", "tpsa", "rotatable_bonds")


def _merge(records: List[CompoundRecord], query: str) -> CompoundRecord:
    base = CompoundRecord(source="merged", source_id=" + ".join(f"{r.source}:{r.source_id}" for r in records), query=query)
    for f in _MERGE_FIELDS:
        for r in records:
            v = getattr(r, f)
            if v not in (None, "", []):
                setattr(base, f, v)
                break
    syn: List[str] = []
    for r in records:
        for s in ([r.name] if r.name else []) + r.synonyms:
            if s and s not in syn:
                syn.append(s)
    base.synonyms = syn[:40]
    for r in records:
        base.extra[r.source] = r.extra
    return base


def resolve(
    query: str,
    sources: Optional[Iterable[str]] = None,
    query_type: str = "auto",
    unichem: bool = True,
    timeout: float = 40.0,
) -> ResolveResult:
    """Look a query up in several databases at once and reconcile the answers.

    Args:
        query: Name, SMILES, InChIKey, or a native identifier (CID, CHEMBL…, CHEBI:…, C…).
        sources: Provider keys; default :data:`DEFAULT_SOURCES`. Use ``["all"]`` for every provider.
        query_type: One of :data:`QUERY_TYPES`; ``auto`` infers from the query's shape.
        unichem: Also fetch cross-database identifiers from UniChem via the consensus InChIKey.
        timeout: Overall wall-clock budget for the parallel fan-out, in seconds.
    """
    if query_type.lower() not in QUERY_TYPES:
        raise ValueError(f"Unknown query type {query_type!r}. Valid: {QUERY_TYPES}")
    keys = list(PROVIDERS) if (sources and "all" in [s.lower() for s in sources]) else [s.lower() for s in (sources or DEFAULT_SOURCES)]
    for k in keys:
        if k not in PROVIDERS:
            raise ValueError(f"Unknown source {k!r}. Available: {sorted(PROVIDERS)}")

    outcomes: List[SourceOutcome] = []
    records: List[CompoundRecord] = []

    deadline = time.time() + timeout

    def _one(k: str):
        t0 = time.time()
        try:
            rec = get_provider(k).lookup(query, query_type, deadline=deadline)
            return k, rec, None, time.time() - t0
        except ProviderError as exc:
            return k, None, str(exc), time.time() - t0
        except Exception as exc:  # defensive: a provider bug must not kill the others
            logger.exception("provider %s failed", k)
            return k, None, f"{type(exc).__name__}: {exc}", time.time() - t0

    with cf.ThreadPoolExecutor(max_workers=len(keys)) as ex:
        futs = {ex.submit(_one, k): k for k in keys}
        done, pending = cf.wait(futs, timeout=timeout)
        for fut in pending:
            outcomes.append(SourceOutcome(source=futs[fut], status="error", error=f"timed out after {timeout:.0f}s", seconds=timeout))
        results = [f.result() for f in done]
    order = {k: i for i, k in enumerate(keys)}
    for k, rec, err, secs in sorted(results, key=lambda x: order[x[0]]):
        if err:
            outcomes.append(SourceOutcome(source=k, status="error", error=err, seconds=round(secs, 2)))
        elif rec is None:
            outcomes.append(SourceOutcome(source=k, status="not found", seconds=round(secs, 2)))
        else:
            outcomes.append(SourceOutcome(source=k, status="found", seconds=round(secs, 2)))
            records.append(rec)
    outcomes.sort(key=lambda o: order[o.source])

    # Second pass: a *name* that one database knows under a different label can be recovered in
    # the others through the structure. Retry every missing source by the consensus InChIKey.
    iks0 = [r.inchikey for r in records if r.inchikey]
    if iks0 and query_type.lower() in ("auto", "name") and not _looks_like_inchikey(query):
        consensus0, n0 = Counter(iks0).most_common(1)[0]
        missing = [o.source for o in outcomes if o.status == "not found"]
        # A source that answered with a *different* structure than the majority most likely
        # matched a synonym or a derivative by text search; re-resolve it by structure.
        outliers = [r.source for r in records if r.inchikey and r.inchikey != consensus0] if n0 >= 2 or len(iks0) == 2 else []
        if outliers:
            records = [r for r in records if r.source not in outliers]
            missing = missing + outliers
        if missing:
            def _retry(k: str):
                t0 = time.time()
                try:
                    return k, get_provider(k).by_inchikey(consensus0), None, time.time() - t0
                except ProviderError as exc:
                    return k, None, str(exc), time.time() - t0
                except Exception as exc:
                    return k, None, f"{type(exc).__name__}: {exc}", time.time() - t0
            with cf.ThreadPoolExecutor(max_workers=len(missing)) as ex:
                futs = {ex.submit(_retry, k): k for k in missing}
                done, _ = cf.wait(futs, timeout=max(5.0, timeout / 2))
                for fut in done:
                    k, rec, err, secs = fut.result()
                    for o in outcomes:
                        if o.source == k:
                            o.seconds = round(o.seconds + secs, 2)
                            if rec is not None:
                                o.status = "found"; o.error = None
                            elif k in outliers:
                                o.status = "not found"; o.error = "text match disagreed with other sources and no structure match was found"
                    if rec is not None:
                        records.append(rec)
            records.sort(key=lambda r: order.get(r.source, 99))

    res = ResolveResult(query=query, query_type=query_type, sources_requested=keys, outcomes=outcomes, records=records)
    if not records:
        return res

    iks = [r.inchikey for r in records if r.inchikey]
    if iks:
        consensus, n = Counter(iks).most_common(1)[0]
        res.consensus_inchikey = consensus
        res.agreement = "unknown" if len(iks) < 2 else ("agree" if n == len(iks) else "disagree")

    xrefs: Dict[str, str] = {}
    for r in records:
        xrefs.setdefault(r.source, r.source_id)
        for k, v in r.cross_refs.items():
            xrefs.setdefault(k, v)
    if unichem and res.consensus_inchikey:
        try:
            for k, v in unichem_xrefs(res.consensus_inchikey).items():
                xrefs.setdefault(k, v)
        except ProviderError as exc:
            logger.debug("UniChem unavailable: %s", exc)
    res.cross_refs = dict(sorted(xrefs.items()))
    res.merged = _merge(records, query)
    res.merged.cross_refs = res.cross_refs
    return res


__all__ = [
    "PROVIDERS", "DEFAULT_SOURCES", "QUERY_TYPES",
    "CompoundRecord", "Provider", "ProviderError", "SourceOutcome", "ResolveResult",
    "get_provider", "resolve", "unichem_xrefs",
]
