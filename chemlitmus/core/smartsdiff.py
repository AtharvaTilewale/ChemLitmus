"""Semantic diff of two SMARTS catalogues.

Two versions of an alert set — or two implementations of the "same" set — are compared on what
they *do* to a reference library, not just on their text:

* patterns are **paired** by name, then by identical SMARTS, then by identical hit set;
* each pair is classified textually (``identical`` / ``rewritten`` / ``no SMARTS``) and
  semantically (``same hits`` / ``broadened`` / ``narrowed`` / ``shifted`` / ``broken``);
* unpaired patterns are ``added`` or ``removed``;
* at catalogue level, the number of library molecules whose **flagged / not flagged verdict
  changes** between the two catalogues is reported — the figure that matters for a screen.

A side can be a pattern file (see :func:`chemlitmus.core.smartsaudit.load_patterns`) or one of
RDKit's built-in filter catalogues written as ``rdkit:<NAME>`` (``rdkit:PAINS``,
``rdkit:PAINS_A``, ``rdkit:BRENK``, ``rdkit:NIH``, ``rdkit:ZINC``, ``rdkit:CHEMBL_Glaxo`` …). RDKit does
not expose the SMARTS text of its catalogues, so those entries are compared on hits alone.

Everything is offline. Hit sets are computed on the chosen molecule preparation (see
``--prep``); a pattern that is dead under one preparation and alive under another is a
*preparation* question (``smartsaudit``), not a catalogue change.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from pydantic import BaseModel, Field

from chemlitmus.core.smartsaudit import PREPARATIONS, _build_library, load_patterns, load_reference_library, prepare_molecule_status

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import FilterCatalog as _fc

    RDLogger.DisableLog("rdApp.*")
    _RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover
    _RDKIT_AVAILABLE = False

RDKIT_CATALOGS: List[str] = [
    "PAINS", "PAINS_A", "PAINS_B", "PAINS_C", "BRENK", "NIH", "ZINC", "CHEMBL", "CHEMBL_BMS", "CHEMBL_Dundee",
    "CHEMBL_Glaxo", "CHEMBL_Inpharmatica", "CHEMBL_LINT", "CHEMBL_MLSMR", "CHEMBL_SureChEMBL", "ALL",
]

TEXT_STATUS = ["identical", "rewritten", "no SMARTS"]
SEMANTIC_STATUS = ["same hits", "broadened", "narrowed", "shifted", "broken", "repaired", "added", "removed"]


class PatternSide(BaseModel):
    name: Optional[str] = None
    smarts: Optional[str] = None
    rule_set: Optional[str] = None
    parses: bool = True
    n_hits: int = 0
    match_error: Optional[str] = None


class PatternDiff(BaseModel):
    """One pattern (or pair of patterns) across the two catalogues."""

    key: str = Field(description="Name when available, else the SMARTS, else a positional id.")
    paired_by: Optional[str] = Field(None, description="'name' | 'smarts' | 'hits' | None (unpaired)")
    a: Optional[PatternSide] = None
    b: Optional[PatternSide] = None
    text_status: Optional[str] = Field(None, description="identical | rewritten | no SMARTS (one side has no text)")
    semantic_status: str = Field(description="same hits | broadened | narrowed | shifted | broken | repaired | added | removed")
    hits_a: int = 0
    hits_b: int = 0
    gained: int = Field(0, description="Molecules hit by B but not A.")
    lost: int = Field(0, description="Molecules hit by A but not B.")
    jaccard: Optional[float] = Field(None, description="|A∩B| / |A∪B| of the hit sets; None when both are empty.")
    examples_gained: List[str] = Field(default_factory=list)
    examples_lost: List[str] = Field(default_factory=list)


class SmartsDiffResult(BaseModel):
    source_a: str
    source_b: str
    library_source: str
    n_molecules: int = Field(description="Reference molecules that reached the requested preparation and were compared.")
    n_unevaluated: int = Field(0, description="Reference molecules excluded because the preparation failed for them.")
    preparation: str
    n_patterns_a: int
    n_patterns_b: int
    n_paired: int
    paired_by: Dict[str, int] = Field(default_factory=dict)
    text_counts: Dict[str, int] = Field(default_factory=dict)
    semantic_counts: Dict[str, int] = Field(default_factory=dict)
    flagged_a: int = Field(0, description="Library molecules flagged by at least one pattern of A.")
    flagged_b: int = 0
    flagged_both: int = 0
    flagged_only_a: int = Field(0, description="Molecules that stop being flagged when moving from A to B.")
    flagged_only_b: int = Field(0, description="Molecules that become flagged when moving from A to B.")
    verdict_changes: int = Field(0, description="flagged_only_a + flagged_only_b: molecules whose pass/fail verdict differs between the catalogues.")
    verdict_change_fraction: float = 0.0
    patterns: List[PatternDiff] = Field(default_factory=list, description="Changed, added and removed patterns first; unchanged ones last.")
    error: Optional[str] = None


# --------------------------------------------------------------------------------------------- loading

def _is_rdkit_source(src: str) -> bool:
    return src.lower().startswith("rdkit:")


def load_side(source: str | Path) -> Tuple[List[PatternSide], List[Optional["Chem.Mol"]], Optional["_fc.FilterCatalog"]]:
    """Patterns of one side: ``(sides, queries, catalog)``.

    For a file, ``queries`` holds compiled SMARTS (``None`` where parsing failed) and ``catalog``
    is ``None``. For ``rdkit:<NAME>``, ``queries`` is empty and ``catalog`` is the FilterCatalog
    whose entries are matched molecule by molecule.
    """
    s = str(source)
    if _is_rdkit_source(s):
        name = s.split(":", 1)[1]
        if name not in RDKIT_CATALOGS:
            raise ValueError(f"Unknown RDKit catalogue {name!r}. Available: {RDKIT_CATALOGS}")
        params = _fc.FilterCatalogParams()
        params.AddCatalog(getattr(_fc.FilterCatalogParams.FilterCatalogs, name))
        cat = _fc.FilterCatalog(params)
        sides = []
        for i in range(cat.GetNumEntries()):
            e = cat.GetEntry(i)
            rs = e.GetProp("FilterSet") if "FilterSet" in list(e.GetPropList()) else name
            sides.append(PatternSide(name=e.GetDescription(), smarts=None, rule_set=rs))
        return sides, [], cat
    pats = load_patterns(s)
    sides, queries = [], []
    for smarts, name, rs in pats:
        q = Chem.MolFromSmarts(smarts)
        sides.append(PatternSide(name=name, smarts=smarts, rule_set=rs, parses=q is not None))
        queries.append(q)
    return sides, queries, None


def _hit_matrix(sides: List[PatternSide], queries, catalog, mols: Sequence["Chem.Mol"], preparation: str) -> np.ndarray:
    n = len(mols)
    M = np.zeros((len(sides), n), dtype=bool)
    evaluated = np.ones(n, dtype=bool)
    if catalog is not None:
        index = {}
        for i in range(catalog.GetNumEntries()):
            index.setdefault(catalog.GetEntry(i).GetDescription(), []).append(i)
        for j, m in enumerate(mols):
            pm, status = prepare_molecule_status(m, preparation)
            if pm is None:
                evaluated[j] = False
                continue
            for e in catalog.GetMatches(pm):
                for i in index.get(e.GetDescription(), []):
                    M[i, j] = True
    else:
        lib = _build_library(mols, preparation)
        evaluated[:] = lib.evaluated
        for i, q in enumerate(queries):
            if q is not None:
                M[i], err = lib.match(q)
                if err:
                    sides[i].parses = False
                    sides[i].match_error = err
    for i, s in enumerate(sides):
        s.n_hits = int(M[i].sum())
    return M, evaluated


# --------------------------------------------------------------------------------------------- pairing

def _pair(sa: List[PatternSide], sb: List[PatternSide], Ma: np.ndarray, Mb: np.ndarray) -> List[Tuple[Optional[int], Optional[int], Optional[str]]]:
    """Return (i_a, i_b, paired_by) triples covering every pattern of both sides."""
    used_a, used_b = set(), set()
    pairs: List[Tuple[Optional[int], Optional[int], Optional[str]]] = []

    def _match(key_fn, label):
        idx_b: Dict[str, List[int]] = {}
        for j, s in enumerate(sb):
            if j in used_b:
                continue
            k = key_fn(s, Mb[j])
            if k is not None:
                idx_b.setdefault(k, []).append(j)
        for i, s in enumerate(sa):
            if i in used_a:
                continue
            k = key_fn(s, Ma[i])
            if k is None or not idx_b.get(k):
                continue
            j = idx_b[k].pop(0)
            used_a.add(i); used_b.add(j)
            pairs.append((i, j, label))

    _match(lambda s, v: s.name.strip().lower() if s.name else None, "name")
    _match(lambda s, v: s.smarts.strip() if s.smarts else None, "smarts")
    # identical non-empty hit sets: renamed / rewritten but behaviourally the same pattern
    _match(lambda s, v: v.tobytes() if v.any() else None, "hits")
    for i in range(len(sa)):
        if i not in used_a:
            pairs.append((i, None, None))
    for j in range(len(sb)):
        if j not in used_b:
            pairs.append((None, j, None))
    return pairs


def _semantic(a: Optional[PatternSide], b: Optional[PatternSide], va, vb) -> str:
    if a is None:
        return "added"
    if b is None:
        return "removed"
    if not a.parses and b.parses:
        return "repaired"
    if a.parses and not b.parses:
        return "broken"
    same = bool(np.array_equal(va, vb))
    if same:
        return "same hits"
    gained = int((vb & ~va).sum()); lost = int((va & ~vb).sum())
    if gained and not lost:
        return "broadened"
    if lost and not gained:
        return "narrowed"
    return "shifted"


# --------------------------------------------------------------------------------------------- public

def diff_smarts(
    source_a: str | Path,
    source_b: str | Path,
    library: Optional[str | Path] = None,
    preparation: str = "implicit-h",
    max_molecules: Optional[int] = None,
    n_examples: int = 3,
) -> SmartsDiffResult:
    """Compare two SMARTS catalogues on a reference library."""
    if not _RDKIT_AVAILABLE:  # pragma: no cover
        raise ImportError("RDKit is required for smartsdiff")
    if preparation not in PREPARATIONS:
        raise ValueError(f"Unknown preparation {preparation!r}. Valid: {PREPARATIONS}")
    mols, lib_src = load_reference_library(library, max_molecules=max_molecules)
    smiles = [Chem.MolToSmiles(m) for m in mols]
    sa, qa, ca = load_side(source_a)
    sb, qb, cb = load_side(source_b)
    Ma, eva = _hit_matrix(sa, qa, ca, mols, preparation)
    Mb, evb = _hit_matrix(sb, qb, cb, mols, preparation)
    evaluated = eva & evb
    n_unevaluated = int((~evaluated).sum())
    # molecules that could not be prepared are excluded from every comparison below
    Ma = Ma[:, evaluated]; Mb = Mb[:, evaluated]
    smiles = [s for s, ok in zip(smiles, evaluated) if ok]
    mols = [m for m, ok in zip(mols, evaluated) if ok]

    diffs: List[PatternDiff] = []
    for i, j, by in _pair(sa, sb, Ma, Mb):
        a = sa[i] if i is not None else None
        b = sb[j] if j is not None else None
        va = Ma[i] if i is not None else np.zeros(len(mols), bool)
        vb = Mb[j] if j is not None else np.zeros(len(mols), bool)
        if a is not None and b is not None:
            if a.smarts is None or b.smarts is None:
                text = "no SMARTS"
            else:
                text = "identical" if a.smarts.strip() == b.smarts.strip() else "rewritten"
        else:
            text = None
        gained_idx = np.flatnonzero(vb & ~va); lost_idx = np.flatnonzero(va & ~vb)
        union = int((va | vb).sum())
        key = (a or b).name or (a or b).smarts or f"#{i if i is not None else j}"
        diffs.append(PatternDiff(
            key=key, paired_by=by, a=a, b=b, text_status=text, semantic_status=_semantic(a, b, va, vb),
            hits_a=int(va.sum()), hits_b=int(vb.sum()), gained=len(gained_idx), lost=len(lost_idx),
            jaccard=(int((va & vb).sum()) / union) if union else None,
            examples_gained=[smiles[k] for k in gained_idx[:n_examples]],
            examples_lost=[smiles[k] for k in lost_idx[:n_examples]],
        ))
    order = {s: k for k, s in enumerate(["broken", "repaired", "shifted", "narrowed", "broadened", "removed", "added", "same hits"])}
    diffs.sort(key=lambda d: (order[d.semantic_status], -(d.gained + d.lost), d.key))

    fa = Ma.any(axis=0) if len(sa) else np.zeros(len(mols), bool)
    fb = Mb.any(axis=0) if len(sb) else np.zeros(len(mols), bool)
    res = SmartsDiffResult(
        source_a=str(source_a), source_b=str(source_b), library_source=lib_src, n_molecules=len(mols), n_unevaluated=n_unevaluated, preparation=preparation,
        n_patterns_a=len(sa), n_patterns_b=len(sb), n_paired=sum(1 for d in diffs if d.paired_by),
        paired_by=dict(Counter(d.paired_by for d in diffs if d.paired_by)),
        text_counts=dict(Counter(d.text_status for d in diffs if d.text_status)),
        semantic_counts=dict(Counter(d.semantic_status for d in diffs)),
        flagged_a=int(fa.sum()), flagged_b=int(fb.sum()), flagged_both=int((fa & fb).sum()),
        flagged_only_a=int((fa & ~fb).sum()), flagged_only_b=int((fb & ~fa).sum()),
        patterns=diffs,
    )
    res.verdict_changes = res.flagged_only_a + res.flagged_only_b
    res.verdict_change_fraction = res.verdict_changes / len(mols) if mols else 0.0
    return res


__all__ = ["diff_smarts", "SmartsDiffResult", "PatternDiff", "PatternSide", "RDKIT_CATALOGS", "load_side"]
