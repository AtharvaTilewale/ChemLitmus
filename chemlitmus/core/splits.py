"""Reproducible, group-aware split generation.

Groups are never broken: every record sharing an identity key (or a scaffold) lands in the same
split. That makes the achieved fractions approximate — a large group cannot be divided — and the
report states the achieved fractions, the group-size distribution, the endpoint balance, and any
constraint that could not be satisfied. Splits are verified with the same leakage implementation
the standalone audit uses.

Strategies: ``random`` (group = the record itself), ``identity`` (group = identity key at the
policy level), ``scaffold`` (group = Murcko scaffold; acyclic molecules are **not** pooled into
one giant group — each is its own group), ``temporal`` and ``source`` (metadata-driven).
"""

from __future__ import annotations

import random
from collections import Counter, defaultdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

from pydantic import BaseModel, Field

from chemlitmus.core.identity import compute_identity
from chemlitmus.core.policy import ChemicalPolicy

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.DisableLog("rdApp.*")
except ImportError:  # pragma: no cover
    Chem = None

STRATEGIES = ["random", "identity", "scaffold", "temporal", "source"]


class SplitAssignment(BaseModel):
    record_id: str
    source_id: Optional[str] = None
    smiles: str
    group_key: str
    split: str


class SplitReport(BaseModel):
    strategy: str
    seed: Optional[int]
    requested_fractions: Dict[str, float]
    achieved_fractions: Dict[str, float]
    n_records: int
    n_excluded: int = Field(0, description="Records without a usable structure (or without the metadata the strategy needs).")
    n_groups: int
    group_sizes: Dict[str, float] = Field(default_factory=dict, description="min / median / max group size and largest_fraction_percent of the dataset.")
    n_by_split: Dict[str, int] = Field(default_factory=dict)
    endpoint_balance: Dict[str, Dict[str, float]] = Field(default_factory=dict, description="split -> {label: fraction} or {mean/min/max} for numeric endpoints.")
    conflicts: List[str] = Field(default_factory=list, description="Constraints that could not all be satisfied, stated rather than silently broken.")
    assignments: List[SplitAssignment] = Field(default_factory=list)
    leakage: Optional[Dict[str, Any]] = Field(None, description="Post-split leakage report, computed with the same implementation as the standalone audit.")
    note: str = ("Groups are indivisible, so achieved fractions differ from requested ones. A scaffold split is an "
                 "evaluation design; it does not guarantee generalisation to deployment data.")


def _scaffold_key(smiles: str, record_id: str) -> str:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return f"invalid:{record_id}"
    try:
        s = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
    except Exception:
        s = ""
    # acyclic molecules have an empty Murcko scaffold; each becomes its own group rather than
    # collapsing every acyclic record into one accidental giant group
    return s if s else f"acyclic:{record_id}"


def _group_keys(records: Sequence[Dict[str, Any]], strategy: str, policy: ChemicalPolicy) -> Tuple[Dict[str, str], List[str]]:
    keys: Dict[str, str] = {}
    excluded: List[str] = []
    for r in records:
        rid, smi = r["record_id"], r.get("smiles")
        if strategy == "random":
            keys[rid] = rid
            continue
        if strategy in ("temporal", "source"):
            field = "date" if strategy == "temporal" else "source"
            v = str(r.get(field, "")).strip()
            if not v:
                excluded.append(rid)
                continue
            keys[rid] = v
            continue
        if not smi:
            excluded.append(rid)
            continue
        if strategy == "scaffold":
            keys[rid] = _scaffold_key(smi, rid)
        else:
            k = compute_identity(smi)
            if not k.is_valid or not k.key(policy.identity_level):
                excluded.append(rid)
                continue
            keys[rid] = k.key(policy.identity_level)
    return keys, excluded


def _endpoint_balance(assignments: List[SplitAssignment], endpoints: Dict[str, str]) -> Dict[str, Dict[str, float]]:
    if not endpoints:
        return {}
    by_split: Dict[str, List[str]] = defaultdict(list)
    for a in assignments:
        v = endpoints.get(a.record_id)
        if v not in (None, ""):
            by_split[a.split].append(str(v))
    out: Dict[str, Dict[str, float]] = {}
    for split, vals in by_split.items():
        try:
            nums = [float(v) for v in vals]
            out[split] = {"n": float(len(nums)), "mean": round(sum(nums) / len(nums), 4), "min": min(nums), "max": max(nums)}
        except ValueError:
            c = Counter(vals)
            out[split] = {lab: round(n / len(vals), 4) for lab, n in c.most_common()}
    return out


def make_splits(
    records: Sequence[Dict[str, Any]],
    fractions: Optional[Dict[str, float]] = None,
    strategy: str = "identity",
    policy: Optional[ChemicalPolicy] = None,
    seed: int = 0,
    endpoint_field: Optional[str] = None,
    verify: bool = True,
) -> SplitReport:
    """Assign records to splits without breaking groups.

    Args:
        records: dicts with ``record_id``, ``smiles`` and, for metadata strategies, ``date`` or ``source``.
        fractions: e.g. ``{"train": 0.8, "test": 0.2}``; must sum to 1.
        strategy: one of :data:`STRATEGIES`.
        seed: makes ``random``/``identity``/``scaffold`` assignments reproducible.
        endpoint_field: key in ``records`` whose distribution is reported per split.
        verify: run the leakage report on the produced splits.
    """
    pol = policy or ChemicalPolicy.preset("parent")
    if strategy not in STRATEGIES:
        raise ValueError(f"Unknown strategy {strategy!r}. Valid: {STRATEGIES}")
    fr = dict(fractions or {"train": 0.8, "test": 0.2})
    if abs(sum(fr.values()) - 1.0) > 1e-6:
        raise ValueError(f"fractions must sum to 1, got {sum(fr.values())}")
    keys, excluded = _group_keys(records, strategy, pol)
    members: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for r in records:
        if r["record_id"] in keys:
            members[keys[r["record_id"]]].append(r)
    n_assignable = sum(len(v) for v in members.values())
    conflicts: List[str] = []
    if excluded:
        conflicts.append(f"{len(excluded)} record(s) excluded: no usable structure or missing {'date' if strategy == 'temporal' else 'source' if strategy == 'source' else 'structure'}")

    names = list(fr)
    order: List[str]
    if strategy in ("temporal", "source"):
        order = sorted(members)                       # chronological / lexical; earliest groups fill the first split
        if strategy == "temporal" and len(names) > 1:
            conflicts.append("temporal split: groups are ordered by date, so the requested fractions are met only to the nearest whole group")
    else:
        rng = random.Random(seed)
        order = sorted(members)
        rng.shuffle(order)
        order.sort(key=lambda k: -len(members[k]))    # place large groups first so they do not overflow the last split
    targets = {n: fr[n] * n_assignable for n in names}
    counts = {n: 0 for n in names}
    assignments: List[SplitAssignment] = []
    for k in order:
        if strategy in ("temporal", "source"):
            # fill splits in order until each reaches its target
            split = names[-1]
            for n in names:
                if counts[n] < targets[n]:
                    split = n
                    break
        else:
            split = min(names, key=lambda n: (counts[n] - targets[n]) / max(targets[n], 1e-9))
        counts[split] += len(members[k])
        for r in members[k]:
            assignments.append(SplitAssignment(record_id=r["record_id"], source_id=r.get("source_id"), smiles=r.get("smiles", ""), group_key=k, split=split))
    sizes = sorted(len(v) for v in members.values())
    largest = sizes[-1] if sizes else 0
    if n_assignable and largest / n_assignable > min(fr.values()):
        conflicts.append(f"largest group holds {largest} records ({largest / n_assignable:.1%}), more than the smallest requested fraction ({min(fr.values()):.1%}); that split is necessarily larger than requested")
    for n in names:
        if n_assignable and abs(counts[n] / n_assignable - fr[n]) > 0.05:
            conflicts.append(f"split {n!r}: achieved {counts[n] / n_assignable:.1%} vs requested {fr[n]:.1%} (indivisible groups)")
    rep = SplitReport(
        strategy=strategy, seed=None if strategy in ("temporal", "source") else seed, requested_fractions=fr,
        achieved_fractions={n: round(counts[n] / n_assignable, 4) if n_assignable else 0.0 for n in names},
        n_records=len(records), n_excluded=len(excluded), n_groups=len(members),
        group_sizes={"min": sizes[0] if sizes else 0, "median": sizes[len(sizes) // 2] if sizes else 0, "max": largest,
                     "largest_fraction_percent": round(100 * largest / n_assignable, 2) if n_assignable else 0},
        n_by_split=counts, conflicts=conflicts, assignments=assignments,
    )
    if endpoint_field:
        rep.endpoint_balance = _endpoint_balance(assignments, {r["record_id"]: r.get(endpoint_field) for r in records})
    if verify:
        from chemlitmus.core.leakage import leakage_report
        splits: Dict[str, List[Tuple[str, str, Optional[str]]]] = defaultdict(list)
        for a in assignments:
            if a.smiles:
                splits[a.split].append((a.record_id, a.smiles, None))
        rep.leakage = leakage_report(splits, pol).model_dump()
    return rep


__all__ = ["make_splits", "SplitReport", "SplitAssignment", "STRATEGIES"]
