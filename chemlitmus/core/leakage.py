"""Chemical leakage between dataset splits.

Overlap is reported as *distinct evidence classes* — exact, parent, tautomer, nostereo,
skeleton — each counted on its own. The classes nest (``exact`` ⊂ ``parent`` ⊂ {tautomer,
nostereo} ⊂ skeleton), so the counts are **not additive**. Scaffold overlap and nearest-neighbour
similarity are *relatedness* diagnostics: they say an evaluation is easy, not that it is invalid.
Formula matches are reported separately and are never called leakage.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

from pydantic import BaseModel, Field

from chemlitmus.core.identity import compute_identity
from chemlitmus.core.policy import ChemicalPolicy

try:
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.DisableLog("rdApp.*")
except ImportError:  # pragma: no cover
    Chem = None

OVERLAP_LEVELS = ["exact", "parent", "tautomer", "nostereo", "skeleton"]


class OverlapClass(BaseModel):
    level: str
    n_eval_records: int = Field(description="Evaluation-split records with a match in the reference split at this level.")
    n_eval_groups: int = Field(description="Distinct keys at this level shared between the two splits.")
    fraction_of_eval: float
    examples: List[Dict[str, Any]] = Field(default_factory=list, description="[{eval_record_id, reference_record_id, key}]")


class NeighbourHit(BaseModel):
    eval_record_id: str
    reference_record_id: str
    similarity: float
    eval_smiles: str
    reference_smiles: str


class PairReport(BaseModel):
    reference_split: str
    evaluation_split: str
    n_reference: int
    n_evaluation: int
    overlap: List[OverlapClass]
    formula_matches: int = Field(0, description="Evaluation records sharing a formula with a reference record. Not identity; reported for completeness.")
    scaffold_overlap_records: int = 0
    scaffold_overlap_fraction: float = 0.0
    n_acyclic_evaluation: int = Field(0, description="Evaluation records without a Murcko scaffold; excluded from scaffold overlap.")
    neighbour_fingerprint: str
    neighbour_threshold: float
    n_related_by_similarity: int = Field(0, description="Evaluation records whose nearest reference neighbour is at or above the threshold.")
    nearest_neighbour_similarity: Dict[str, float] = Field(default_factory=dict, description="Summary of the nearest-neighbour Tanimoto over evaluation records: min/median/max.")
    neighbours: List[NeighbourHit] = Field(default_factory=list, description="Per evaluation record: its nearest reference neighbour.")
    temporal_violations: int = Field(0, description="Evaluation records dated no later than the latest reference record (when dates are available).")
    temporal_note: Optional[str] = None


class LeakageReport(BaseModel):
    identity_level: str
    splits: List[str]
    n_by_split: Dict[str, int]
    n_unassigned: int = Field(0, description="Records without a split label (excluded).")
    within_split_duplicates: Dict[str, int] = Field(default_factory=dict, description="Per split: records collapsed at the policy identity level.")
    pairs: List[PairReport]
    note: str = ("Overlap classes nest and are not additive. Scaffold and similarity relatedness are evaluation-design "
                 "diagnostics, not evidence of an invalid experiment. Formula matches are not identity.")
    record_issues: Dict[str, List[Dict[str, Any]]] = Field(default_factory=dict, description="record_id -> issue dicts (code + evidence) for the audit to attach.")


# ------------------------------------------------------------------------------------ core

def _fp(smiles: str, policy: ChemicalPolicy):
    from chemlitmus.core.cheminfo import _compute_single_fp
    m = Chem.MolFromSmiles(smiles)
    return _compute_single_fp(m, policy.fingerprint, policy.fingerprint_bits)[0] if m is not None else None


def _scaffold(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return "acyclic"
    try:
        s = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
    except Exception:
        s = ""
    return s or "acyclic"


class _Item:
    __slots__ = ("record_id", "smiles", "keys", "scaffold", "fp", "date")

    def __init__(self, record_id, smiles, keys, scaffold, fp, date=None):
        self.record_id, self.smiles, self.keys, self.scaffold, self.fp, self.date = record_id, smiles, keys, scaffold, fp, date


def _items(records: Sequence[Tuple[str, str, Optional[str]]], policy: ChemicalPolicy) -> List[_Item]:
    out = []
    for rid, smi, date in records:
        k = compute_identity(smi)
        if not k.is_valid:
            continue
        out.append(_Item(rid, smi, k, _scaffold(smi), _fp(smi, policy), date))
    return out


def _pair(ref: List[_Item], ev: List[_Item], ref_name: str, ev_name: str, policy: ChemicalPolicy, issues: Dict[str, List[Dict[str, Any]]], n_examples: int = 10) -> PairReport:
    overlap: List[OverlapClass] = []
    for lvl in OVERLAP_LEVELS:
        ref_keys: Dict[str, str] = {}
        for r in ref:
            k = r.keys.key(lvl)
            if k and k not in ref_keys:
                ref_keys[k] = r.record_id
        hits = [(e, ref_keys[e.keys.key(lvl)]) for e in ev if e.keys.key(lvl) in ref_keys]
        overlap.append(OverlapClass(level=lvl, n_eval_records=len(hits), n_eval_groups=len({e.keys.key(lvl) for e, _ in hits}),
                                    fraction_of_eval=round(len(hits) / len(ev), 4) if ev else 0.0,
                                    examples=[{"eval_record_id": e.record_id, "reference_record_id": rr, "key": e.keys.key(lvl)} for e, rr in hits[:n_examples]]))
        if lvl == policy.identity_level:
            for e, rr in hits:
                issues.setdefault(e.record_id, []).append({"code": "SPLIT_OVERLAP", "level": lvl, "reference_split": ref_name, "reference_record_id": rr})
    ref_formulas = {r.keys.formula for r in ref if r.keys.formula}
    formula_matches = sum(1 for e in ev if e.keys.formula in ref_formulas)
    ref_scaf = {r.scaffold for r in ref if r.scaffold != "acyclic"}
    cyclic_ev = [e for e in ev if e.scaffold != "acyclic"]
    scaf_hits = [e for e in cyclic_ev if e.scaffold in ref_scaf]
    # nearest neighbour
    ref_fps = [r.fp for r in ref if r.fp is not None]
    ref_ids = [r.record_id for r in ref if r.fp is not None]
    neighbours: List[NeighbourHit] = []
    for e in ev:
        if e.fp is None or not ref_fps:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(e.fp, ref_fps)
        j = max(range(len(sims)), key=sims.__getitem__)
        neighbours.append(NeighbourHit(eval_record_id=e.record_id, reference_record_id=ref_ids[j], similarity=round(sims[j], 4),
                                       eval_smiles=e.smiles, reference_smiles=ref[[r.record_id for r in ref].index(ref_ids[j])].smiles))
    sims_all = sorted(n.similarity for n in neighbours)
    nn_summary = {}
    if sims_all:
        nn_summary = {"min": sims_all[0], "median": sims_all[len(sims_all) // 2], "max": sims_all[-1]}
    related = [n for n in neighbours if n.similarity >= policy.similarity_threshold]
    exact_ids = {ex["eval_record_id"] for oc in overlap for ex in oc.examples}
    for n in related:
        if n.eval_record_id not in [i for i in issues if any(d["code"] == "SPLIT_OVERLAP" for d in issues[i])]:
            issues.setdefault(n.eval_record_id, []).append({"code": "SPLIT_RELATED", "kind": "similarity", "similarity": n.similarity, "reference_record_id": n.reference_record_id, "reference_split": ref_name})
    for e in scaf_hits:
        if not any(d["code"] in ("SPLIT_OVERLAP", "SPLIT_RELATED") for d in issues.get(e.record_id, [])):
            issues.setdefault(e.record_id, []).append({"code": "SPLIT_RELATED", "kind": "scaffold", "scaffold": e.scaffold, "reference_split": ref_name})
    del exact_ids
    # temporal
    t_viol, t_note = 0, None
    ref_dates = [r.date for r in ref if r.date]
    ev_dates = [(e, e.date) for e in ev if e.date]
    if ref_dates and ev_dates:
        latest_ref = max(ref_dates)
        for e, d in ev_dates:
            if d <= latest_ref:
                t_viol += 1
                issues.setdefault(e.record_id, []).append({"code": "TEMPORAL_ORDER", "date": d, "latest_reference_date": latest_ref, "reference_split": ref_name})
        missing = len(ref) - len(ref_dates) + len(ev) - len(ev_dates)
        t_note = f"dates compared as strings (ISO-like ordering assumed); {missing} records without a date were not compared"
    else:
        t_note = "no dates available in both splits; temporal check not run"
    return PairReport(
        reference_split=ref_name, evaluation_split=ev_name, n_reference=len(ref), n_evaluation=len(ev), overlap=overlap,
        formula_matches=formula_matches, scaffold_overlap_records=len(scaf_hits), scaffold_overlap_fraction=round(len(scaf_hits) / len(cyclic_ev), 4) if cyclic_ev else 0.0,
        n_acyclic_evaluation=len(ev) - len(cyclic_ev), neighbour_fingerprint=f"{policy.fingerprint}/{policy.fingerprint_bits}", neighbour_threshold=policy.similarity_threshold,
        n_related_by_similarity=len(related), nearest_neighbour_similarity=nn_summary, neighbours=neighbours, temporal_violations=t_viol, temporal_note=t_note,
    )


def leakage_report(
    splits: Dict[str, Sequence[Tuple[str, str, Optional[str]]]],
    policy: Optional[ChemicalPolicy] = None,
    reference: Optional[str] = None,
) -> LeakageReport:
    """Leakage between named splits.

    Args:
        splits: ``{split_name: [(record_id, smiles, date_or_None), ...]}``.
        policy: Identity level, fingerprint and similarity threshold come from here.
        reference: The split treated as training; default the largest (ties: first). Every other
            split is evaluated against it, and against each other pairwise in the report order.
    """
    pol = policy or ChemicalPolicy.preset("parent")
    items = {name: _items(recs, pol) for name, recs in splits.items()}
    names = list(items)
    if not names:
        return LeakageReport(identity_level=pol.identity_level, splits=[], n_by_split={}, pairs=[])
    ref_name = reference or max(names, key=lambda n: (len(items[n]), -names.index(n)))
    issues: Dict[str, List[Dict[str, Any]]] = {}
    pairs = [_pair(items[ref_name], items[n], ref_name, n, pol, issues) for n in names if n != ref_name]
    others = [n for n in names if n != ref_name]
    for i in range(len(others)):
        for j in range(i + 1, len(others)):
            pairs.append(_pair(items[others[i]], items[others[j]], others[i], others[j], pol, issues))
    within = {}
    for n, its in items.items():
        keys = [it.keys.key(pol.identity_level) for it in its if it.keys.key(pol.identity_level)]
        within[n] = len(keys) - len(set(keys))
    return LeakageReport(identity_level=pol.identity_level, splits=names, n_by_split={n: len(items[n]) for n in names},
                         within_split_duplicates=within, pairs=pairs, record_issues=issues)


def leakage_from_annotations(anns, policy: ChemicalPolicy, groups=None, date_field: Optional[str] = None) -> LeakageReport:
    """Adapter used by the dataset audit (annotations already carry split labels and fields)."""
    splits: Dict[str, List[Tuple[str, str, Optional[str]]]] = defaultdict(list)
    unassigned = 0
    for a in anns:
        if not a.split:
            unassigned += 1
            continue
        date = str(a.fields.get(date_field)).strip() if date_field and a.fields.get(date_field) not in (None, "") else None
        splits[a.split].append((a.record_id, a.parsed_smiles, date))
    rep = leakage_report(splits, policy)
    rep.n_unassigned = unassigned
    return rep


__all__ = ["leakage_report", "leakage_from_annotations", "LeakageReport", "PairReport", "OverlapClass", "NeighbourHit", "OVERLAP_LEVELS"]
