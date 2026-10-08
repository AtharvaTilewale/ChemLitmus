"""Name-to-structure concordance across databases.

Given a list of queries (typically compound names), resolve each across several databases and
tabulate how often the databases return the *same structure* for the *same name* — and, when
they do not, at which identity level they part ways (salt form, tautomer, stereochemistry, or a
different compound altogether). Also counts how often a database's text search landed on a
different structure than the majority and was corrected through the consensus InChIKey.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Callable, Dict, Iterable, List, Optional

from pydantic import BaseModel, Field

from chemlitmus.providers import Correction, ResolveResult, resolve


class ConcordanceRow(BaseModel):
    query: str
    n_found: int
    agreement: str
    agreement_level: Optional[str] = None
    disagreement: Optional[str] = None
    consensus_inchikey: Optional[str] = None
    ids: Dict[str, str] = Field(default_factory=dict, description="source -> native identifier")
    names: Dict[str, str] = Field(default_factory=dict, description="source -> preferred name")
    inchikeys: Dict[str, str] = Field(default_factory=dict)
    corrections: List[Correction] = Field(default_factory=list)
    errors: Dict[str, str] = Field(default_factory=dict, description="source -> error or 'not found'")


class ConcordanceReport(BaseModel):
    sources: List[str]
    n_queries: int
    n_found_any: int = 0
    n_found_all: int = 0
    agreement_counts: Dict[str, int] = Field(default_factory=dict, description="agree / disagree / unknown / not found")
    level_counts: Dict[str, int] = Field(default_factory=dict, description="Strictest identity level shared by all returned structures, for queries with at least two usable structures ('none' is reported as 'different compounds').")
    disagreement_classes: Dict[str, int] = Field(default_factory=dict, description="What differed, for queries with more than one structure returned.")
    found_by_source: Dict[str, int] = Field(default_factory=dict)
    errors_by_source: Dict[str, int] = Field(default_factory=dict)
    corrections_by_source: Dict[str, int] = Field(default_factory=dict, description="Text hits that disagreed with the majority and were re-resolved by structure.")
    correction_classes_by_source: Dict[str, Dict[str, int]] = Field(default_factory=dict)
    rows: List[ConcordanceRow] = Field(default_factory=list)


def row_from_result(res: ResolveResult) -> ConcordanceRow:
    row = ConcordanceRow(query=res.query, n_found=len(res.records), agreement=res.agreement if res.records else "not found",
                         agreement_level=res.agreement_level, disagreement=res.disagreement,
                         consensus_inchikey=res.consensus_inchikey, corrections=res.corrections)
    for r in res.records:
        row.ids[r.source] = r.source_id
        if r.name:
            row.names[r.source] = r.name
        if r.inchikey:
            row.inchikeys[r.source] = r.inchikey
    for o in res.outcomes:
        if o.status != "found":
            row.errors[o.source] = o.error or o.status
    return row


def concordance(
    queries: Iterable[str],
    sources: Optional[List[str]] = None,
    query_type: str = "auto",
    timeout: float = 40.0,
    use_cache: bool = True,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> ConcordanceReport:
    """Resolve every query and tabulate cross-database agreement."""
    qs = [q for q in queries if q and q.strip()]
    rows: List[ConcordanceRow] = []
    srcs: List[str] = []
    for n, q in enumerate(qs, 1):
        res = resolve(q, sources=sources, query_type=query_type, unichem=False, timeout=timeout, use_cache=use_cache)
        srcs = res.sources_requested
        rows.append(row_from_result(res))
        if progress_callback:
            progress_callback(n, len(qs))
    return summarise(rows, srcs)


def summarise(rows: List[ConcordanceRow], sources: List[str]) -> ConcordanceReport:
    rep = ConcordanceReport(sources=sources, n_queries=len(rows), rows=rows)
    rep.n_found_any = sum(1 for r in rows if r.n_found)
    rep.n_found_all = sum(1 for r in rows if r.n_found == len(sources))
    rep.agreement_counts = dict(Counter(r.agreement for r in rows))
    rep.level_counts = dict(Counter(("different compounds" if r.agreement_level == "none" else r.agreement_level)
                                    for r in rows if r.agreement_level is not None))
    rep.disagreement_classes = dict(Counter(r.disagreement for r in rows if r.disagreement))
    found, errs, corr = Counter(), Counter(), Counter()
    classes: Dict[str, Counter] = defaultdict(Counter)
    for r in rows:
        for s in r.ids:
            found[s] += 1
        for s, e in r.errors.items():
            if e != "not found":
                errs[s] += 1
        for c in r.corrections:
            corr[c.source] += 1
            classes[c.source][c.difference] += 1
    rep.found_by_source = {s: found.get(s, 0) for s in sources}
    rep.errors_by_source = {s: errs.get(s, 0) for s in sources}
    rep.corrections_by_source = {s: corr.get(s, 0) for s in sources}
    rep.correction_classes_by_source = {s: dict(c) for s, c in classes.items()}
    return rep


__all__ = ["ConcordanceRow", "ConcordanceReport", "concordance", "summarise", "row_from_result"]
