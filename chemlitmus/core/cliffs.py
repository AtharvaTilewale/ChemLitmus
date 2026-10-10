"""Activity cliffs and label outliers among structurally near-identical compounds.

``conflicts`` compares records of the *same* compound. This module compares *different*
compounds that are nearly the same — matched molecular pairs (one site changed, Hussain & Rea
2010) and high-similarity neighbours — and asks whether their labels agree.

A large label difference across a small structural change is an **activity cliff**: real
structure–activity information when both measurements are sound, a data error when one is not.
The module does not decide which. It reports every pair with its evidence (the structural
relationship, both measurements and the smallest difference the measurements *prove*), and
separately flags **label outliers** — compounds that disagree with several neighbours which
agree with each other — because that is the pattern a wrong value produces most often.

Semantics:

* compounds are compared at the policy identity level; records of one compound are first reduced
  to one value per compound, and compounds whose own records disagree are excluded here (they
  belong to ``conflicts``);
* censored measurements are bounds: a pair is a *cliff* only when the bounds prove a difference at
  or above the threshold, *consistent* only when they prove a difference below it, and
  *undetermined* otherwise — never averaged, never assumed;
* pairs are compared only within the same endpoint context (target / assay columns);
* the molar scale is log10 (pX), so a threshold of 1.0 means ten-fold; non-molar quantities are
  compared in their own units.
"""

from __future__ import annotations

import math
import statistics
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from pydantic import BaseModel, Field
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMMPA

from chemlitmus.core.identity import compute_identity
from chemlitmus.core.labels import CENSOR_RE, parse_measurement
from chemlitmus.core.policy import ChemicalPolicy

VERDICTS = ["cliff", "consistent", "undetermined"]
DEFAULT_MAX_R_ATOMS = 13
DEFAULT_SIMILARITY = 0.9
H_FRAGMENT = "[H][*:1]"


# --------------------------------------------------------------------------- models
class CompoundValue(BaseModel):
    """One compound (identity group) with the value its records support."""

    compound_id: str
    key: str
    smiles: str = Field(description="Representative structure (parent form) used for pairing.")
    record_ids: List[str]
    source_ids: List[str] = Field(default_factory=list, description="User-supplied identifiers of the records, when an id column exists.")
    label: Optional[str] = Field(None, description="Classification label.")
    lo: Optional[float] = Field(None, description="Lower bound on the comparison scale; None = unbounded.")
    hi: Optional[float] = Field(None, description="Upper bound on the comparison scale; None = unbounded.")
    display: str = Field("", description="Human-readable value, e.g. 'pIC50 7.2' or '< 10 uM'.")
    n_exact: int = 0
    n_censored: int = 0

    @property
    def exact(self) -> Optional[float]:
        return self.lo if self.lo is not None and self.lo == self.hi else None


class CliffPair(BaseModel):
    pair_id: str
    context: Dict[str, str] = Field(default_factory=dict)
    compound_a: str
    compound_b: str
    record_ids_a: List[str]
    record_ids_b: List[str]
    source_ids_a: List[str] = Field(default_factory=list)
    source_ids_b: List[str] = Field(default_factory=list)
    smiles_a: str
    smiles_b: str
    relationship: str = Field(description="mmp | similarity | mmp+similarity")
    similarity: Optional[float] = None
    core: Optional[str] = Field(None, description="Shared MMP core with the attachment point as [*:1].")
    transformation: Optional[str] = Field(None, description="R_a>>R_b, the fragments exchanged between a and b.")
    n_changed_atoms: Optional[int] = Field(None, description="Heavy atoms in the two exchanged fragments together.")
    value_a: str
    value_b: str
    min_difference: Optional[float] = Field(None, description="Smallest difference the two measurements prove (0 when their bounds overlap).")
    max_difference: Optional[float] = Field(None, description="Largest difference compatible with the measurements; None when a bound is open.")
    signed_difference: Optional[float] = Field(None, description="value_b - value_a when both are exact.")
    verdict: str = Field(description="cliff | consistent | undetermined")


class TransformationSummary(BaseModel):
    transformation: str = Field(description="Canonical orientation R1>>R2 (R1 sorts first).")
    n_pairs: int
    n_cliffs: int
    n_consistent: int
    n_undetermined: int
    mean_signed_difference: Optional[float] = Field(None, description="Mean of value(R2 compound) - value(R1 compound) over pairs with exact values.")
    n_exact_pairs: int = 0


class LabelOutlier(BaseModel):
    compound_id: str
    record_ids: List[str]
    source_ids: List[str] = Field(default_factory=list)
    smiles: str
    value: str
    context: Dict[str, str] = Field(default_factory=dict)
    n_neighbours: int = Field(description="Near neighbours, every one of which is a proven cliff against this compound.")
    n_agreeing: int = Field(description="Largest set of those neighbours that are proven consistent with one another.")
    neighbour_ids: List[str]
    neighbour_values: List[str]
    neighbour_relationships: List[str]
    note: str = "Every near neighbour disagrees with this compound, none of them disagree with each other, and several agree outright. A wrong value produces this pattern; so does a genuine cliff. Review the source record before changing anything."


class CliffReport(BaseModel):
    endpoint_field: str
    context_fields: List[str]
    identity_level: str
    kind: str = Field(description="classification | quantitative")
    threshold: float
    scale: Optional[str] = Field(None, description="'log10 molar' or 'absolute (<units>)' for quantitative endpoints.")
    use_mmp: bool
    use_similarity: bool
    similarity_threshold: float
    fingerprint: str
    max_r_atoms: int
    n_records: int = Field(description="Records with a parseable structure and a non-empty endpoint value.")
    n_no_value: int = Field(0, description="Parseable records whose endpoint field is empty; not compared.")
    n_compounds: int = Field(description="Distinct compounds at the identity level with a usable value.")
    n_internal_conflict_excluded: int = Field(0, description="Compounds whose own records disagree beyond the threshold; see `conflicts`.")
    n_unparseable: int = 0
    n_missing_units: int = 0
    n_pairs: int
    n_pairs_by_relationship: Dict[str, int]
    n_cliffs: int
    n_consistent: int
    n_undetermined: int
    cliff_fraction: Optional[float] = Field(None, description="n_cliffs / (n_cliffs + n_consistent); undetermined pairs are not in the denominator.")
    n_compounds_in_cliffs: int
    n_outliers: int
    pairs: List[CliffPair] = Field(default_factory=list)
    transformations: List[TransformationSummary] = Field(default_factory=list)
    outliers: List[LabelOutlier] = Field(default_factory=list)
    record_issues: Dict[str, List[Dict[str, Any]]] = Field(default_factory=dict)
    note: str = ("A cliff is a proven label difference across a small structural change: SAR when both values are sound, an error when one is not. "
                 "Censored values are bounds; pairs they cannot decide are 'undetermined', not cliffs. Records of the same compound are compared by `conflicts`, not here.")


# --------------------------------------------------------------------------- intervals
def _interval(m, scale: str) -> Tuple[Optional[float], Optional[float], bool]:
    """Interval on the comparison scale implied by one measurement. Returns (lo, hi, usable)."""
    if m.value is None or not m.comparable:
        return None, None, False
    if scale == "log10 molar":
        if m.value_nm is None or m.value_nm <= 0:
            return None, None, False
        p = -math.log10(m.value_nm * 1e-9)
        if m.relation in ("=", "~"):
            return p, p, True
        if m.relation in ("<", "<="):          # concentration below x  ->  pX above p(x)
            return p, None, True
        if m.relation in (">", ">="):
            return None, p, True
        if m.relation == "range" and m.upper is not None and m.upper > 0:
            factor = m.value_nm / m.value
            p_hi = -math.log10(m.upper * factor * 1e-9)
            return min(p, p_hi), max(p, p_hi), True
        return None, None, False
    v = m.value
    if m.relation in ("=", "~"):
        return v, v, True
    if m.relation in ("<", "<="):
        return None, v, True
    if m.relation in (">", ">="):
        return v, None, True
    if m.relation == "range" and m.upper is not None:
        return min(v, m.upper), max(v, m.upper), True
    return None, None, False


def _gap(a: CompoundValue, b: CompoundValue) -> Tuple[float, Optional[float]]:
    """(minimum proven difference, maximum compatible difference or None if unbounded)."""
    terms = []
    if a.lo is not None and b.hi is not None:
        terms.append(a.lo - b.hi)
    if b.lo is not None and a.hi is not None:
        terms.append(b.lo - a.hi)
    mn = max([0.0] + terms)
    if None in (a.lo, a.hi, b.lo, b.hi):
        return mn, None
    return mn, max(a.hi - b.lo, b.hi - a.lo)


def _verdict(a: CompoundValue, b: CompoundValue, kind: str, threshold: float) -> Tuple[str, Optional[float], Optional[float], Optional[float]]:
    if kind == "classification":
        return ("cliff" if a.label != b.label else "consistent"), None, None, None
    mn, mx = _gap(a, b)
    signed = (b.exact - a.exact) if (a.exact is not None and b.exact is not None) else None
    if mn >= threshold:
        return "cliff", mn, mx, signed
    if mx is not None and mx < threshold:
        return "consistent", mn, mx, signed
    return "undetermined", mn, mx, signed


# --------------------------------------------------------------------------- compounds
def _fmt(v: float) -> str:
    return f"{v:.4g}"


def _compound_values(members: Sequence[Dict[str, Any]], endpoint_field: str, units_field: Optional[str], relation_field: Optional[str],
                     kind: str, threshold: float, counters: Dict[str, int]) -> Tuple[Optional[str], Optional[float], Optional[float], str, int, int, Optional[str]]:
    """Reduce one compound's records to (label, lo, hi, display, n_exact, n_censored, scale) or mark it conflicting/unusable.

    Returns label=None and lo=hi=None with display '' when nothing usable; sets counters['conflict'] when
    the records disagree.
    """
    if kind == "classification":
        labels = {str(r[endpoint_field]).strip() for r in members}
        if len(labels) > 1:
            counters["conflict"] += 1
            return None, None, None, "", 0, 0, None
        lab = labels.pop()
        return lab, None, None, lab, 0, 0, None
    ms = []
    for r in members:
        m = parse_measurement(r[endpoint_field], r.get(units_field) if units_field else None, r.get(relation_field) if relation_field else None, unitless=units_field is None)
        if m.value is None:
            counters["unparseable"] += 1
            continue
        if units_field and not m.units:
            counters["missing_units"] += 1
        ms.append(m)
    usable = [m for m in ms if m.comparable]
    if not usable:
        return None, None, None, "", 0, 0, None
    molar = [m for m in usable if m.log_value is not None or (m.units and m.value_nm is not None)]
    if len(molar) == len(usable):
        scale = "log10 molar"
    else:
        units = {m.units.lower() for m in usable if m.units}
        if len(units) == 1:
            scale = f"absolute ({units.pop()})"
        elif not units:
            scale = "absolute (unitless)"
        else:
            return None, None, None, "", 0, 0, None
    ivs = [(m, *_interval(m, scale)) for m in usable]
    ivs = [(m, lo, hi) for m, lo, hi, ok in ivs if ok]
    if not ivs:
        return None, None, None, "", 0, 0, scale
    exact = [(m, lo) for m, lo, hi in ivs if lo is not None and lo == hi]
    cens = [(m, lo, hi) for m, lo, hi in ivs if not (lo is not None and lo == hi)]
    unit_label = "p" if scale == "log10 molar" else ""
    if exact:
        pts = [lo for _, lo in exact]
        if max(pts) - min(pts) > threshold:
            counters["conflict"] += 1
            return None, None, None, "", len(exact), len(cens), scale
        v = statistics.median(pts)
        disp = f"{unit_label}{_fmt(v)}" + (f" (median of {len(pts)})" if len(pts) > 1 else "")
        if scale != "log10 molar":
            disp += f" {exact[0][0].units or ''}".rstrip()
        return None, v, v, disp, len(exact), len(cens), scale
    lo = max([l for _, l, _ in cens if l is not None], default=None)
    hi = min([h for _, _, h in cens if h is not None], default=None)
    if lo is not None and hi is not None and lo > hi:
        counters["conflict"] += 1
        return None, None, None, "", 0, len(cens), scale
    if lo is None and hi is None:
        return None, None, None, "", 0, len(cens), scale
    if lo is not None and hi is not None:
        disp = f"{unit_label}{_fmt(lo)}–{_fmt(hi)}"
    elif lo is not None:
        disp = f"{unit_label}≥ {_fmt(lo)}"
    else:
        disp = f"{unit_label}≤ {_fmt(hi)}"
    if len(cens) == 1:
        m = cens[0][0]
        raw = m.raw_value if CENSOR_RE.match(m.raw_value) and CENSOR_RE.match(m.raw_value).group(1) else f"{m.relation} {m.raw_value}"
        disp += f" [{raw} {m.units or ''}]".rstrip()
    elif len(cens) > 1:
        disp += f" [{len(cens)} censored]"
    return None, lo, hi, disp, 0, len(cens), scale


# --------------------------------------------------------------------------- pairing
def mmp_fragments(mol: Chem.Mol, max_r_atoms: int = DEFAULT_MAX_R_ATOMS, include_hydrogen: bool = True) -> List[Tuple[str, str, int]]:
    """Single-cut matched-molecular-pair fragments: ``(core, r_group, r_heavy_atoms)``.

    The varied part is at most ``max_r_atoms`` heavy atoms and never larger than the core.
    With ``include_hydrogen`` each H-bearing atom also yields the whole molecule as a core with
    an ``[H]`` R-group, so H -> substituent pairs are found.
    """
    out = set()
    try:
        frags = rdMMPA.FragmentMol(mol, minCuts=1, maxCuts=1, maxCutBonds=40, resultsAsMols=False)
    except Exception:
        frags = ()
    for _core, chains in frags:
        parts = chains.split(".")
        if len(parts) != 2:
            continue
        mols = [Chem.MolFromSmiles(p) for p in parts]
        if any(m is None for m in mols):
            continue
        heavy = [m.GetNumHeavyAtoms() for m in mols]
        for ci, ri in ((0, 1), (1, 0)):
            if heavy[ri] <= max_r_atoms and heavy[ri] <= heavy[ci]:
                out.add((Chem.MolToSmiles(mols[ci]), Chem.MolToSmiles(mols[ri]), heavy[ri]))
    if include_hydrogen:
        # Replacing one H by the attachment point changes nothing but that atom's H count, so
        # re-sanitising is only needed when stereo descriptors could be re-perceived.
        has_stereo = any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or \
            any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds())
        for atom in mol.GetAtoms():
            if atom.GetTotalNumHs() == 0:
                continue
            rw = Chem.RWMol(mol)
            idx = rw.AddAtom(Chem.Atom(0))
            rw.GetAtomWithIdx(idx).SetAtomMapNum(1)
            rw.AddBond(atom.GetIdx(), idx, Chem.BondType.SINGLE)
            a = rw.GetAtomWithIdx(atom.GetIdx())
            if a.GetNumExplicitHs() > 0:
                a.SetNumExplicitHs(a.GetNumExplicitHs() - 1)
            try:
                if has_stereo:
                    Chem.SanitizeMol(rw)
                else:
                    a.UpdatePropertyCache(strict=False)
                out.add((Chem.MolToSmiles(rw), H_FRAGMENT, 0))
            except Exception:
                continue
    return sorted(out)


def _fp(mol, policy: ChemicalPolicy):
    from chemlitmus.core.cheminfo import _compute_single_fp
    return _compute_single_fp(mol, policy.fingerprint, policy.fingerprint_bits)[0]


def _find_pairs(comps: List[CompoundValue], policy: ChemicalPolicy, use_mmp: bool, use_similarity: bool, sim_threshold: float, max_r_atoms: int):
    """Return {(i, j): dict(relationship, similarity, core, transformation, n_changed_atoms)} with i < j."""
    mols = [Chem.MolFromSmiles(c.smiles) for c in comps]
    pairs: Dict[Tuple[int, int], Dict[str, Any]] = {}
    if use_mmp:
        index: Dict[str, Dict[int, set]] = defaultdict(lambda: defaultdict(set))
        for i, m in enumerate(mols):
            if m is None:
                continue
            for core, r, rh in mmp_fragments(m, max_r_atoms):
                index[core][i].add((r, rh))
        for core, members in index.items():
            if len(members) < 2:
                continue
            ids = sorted(members)
            core_heavy = Chem.MolFromSmiles(core).GetNumHeavyAtoms()
            for x in range(len(ids)):
                for y in range(x + 1, len(ids)):
                    i, j = ids[x], ids[y]
                    best = None
                    for ra, ha in members[i]:
                        for rb, hb in members[j]:
                            if ra == rb:
                                continue
                            cand = (ha + hb, ra, rb)
                            if best is None or cand < best:
                                best = cand
                    if best is None:
                        continue
                    prev = pairs.get((i, j))
                    # keep the largest shared core (smallest change)
                    if prev is None or prev.get("core") is None or best[0] < prev["n_changed_atoms"] or (best[0] == prev["n_changed_atoms"] and core_heavy > prev["_core_heavy"]):
                        pairs[(i, j)] = {"relationship": "mmp", "similarity": prev.get("similarity") if prev else None, "core": core,
                                         "transformation": f"{best[1]}>>{best[2]}", "n_changed_atoms": best[0], "_core_heavy": core_heavy}
    if use_similarity:
        fps = [(_fp(m, policy) if m is not None else None) for m in mols]
        valid = [i for i, f in enumerate(fps) if f is not None]
        for a in range(len(valid)):
            i = valid[a]
            rest = valid[a + 1:]
            if not rest:
                continue
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], [fps[j] for j in rest])
            for j, s in zip(rest, sims):
                if s >= sim_threshold:
                    d = pairs.get((i, j))
                    if d is None:
                        pairs[(i, j)] = {"relationship": "similarity", "similarity": round(s, 4), "core": None, "transformation": None, "n_changed_atoms": None}
                    else:
                        d["relationship"] = "mmp+similarity"
                        d["similarity"] = round(s, 4)
        # similarity for mmp-only pairs too, for information
        for (i, j), d in pairs.items():
            if d["similarity"] is None and fps[i] is not None and fps[j] is not None:
                d["similarity"] = round(DataStructs.TanimotoSimilarity(fps[i], fps[j]), 4)
    for d in pairs.values():
        d.pop("_core_heavy", None)
    return pairs


def _largest_consistent_subset(js: List[int], verdicts: Dict[Tuple[int, int], str]) -> int:
    """Size of the largest subset of ``js`` whose members are pairwise 'consistent' (exact for <= 12 members, greedy beyond)."""
    ok = lambda a, b: verdicts.get((a, b)) == "consistent"  # noqa: E731
    if len(js) <= 12:
        best = 0
        n = len(js)
        for mask in range(1, 1 << n):
            sub = [js[k] for k in range(n) if mask >> k & 1]
            if len(sub) <= best:
                continue
            if all(ok(sub[x], sub[y]) for x in range(len(sub)) for y in range(x + 1, len(sub))):
                best = len(sub)
        return best
    chosen: List[int] = []
    for j in js:
        if all(ok(j, k) for k in chosen):
            chosen.append(j)
    return len(chosen)


# --------------------------------------------------------------------------- main
def activity_cliffs(
    records: Sequence[Dict[str, Any]],
    policy: Optional[ChemicalPolicy] = None,
    endpoint_field: str = "endpoint",
    units_field: Optional[str] = None,
    relation_field: Optional[str] = None,
    context_fields: Sequence[str] = (),
    threshold: float = 1.0,
    kind: str = "auto",
    use_mmp: bool = True,
    use_similarity: bool = True,
    similarity_threshold: float = DEFAULT_SIMILARITY,
    max_r_atoms: int = DEFAULT_MAX_R_ATOMS,
    min_outlier_neighbours: int = 2,
    keep_pairs: str = "all",
) -> CliffReport:
    """Find activity cliffs and label outliers.

    Args:
        records: dicts with ``record_id``, ``smiles`` (a parseable structure), optional ``source_id``,
            and the raw field values. Identity keys are computed here at the policy level.
        threshold: label difference that counts as a cliff — log10 units for molar quantities
            (1.0 = ten-fold), absolute otherwise. The same value bounds the disagreement tolerated
            among records of one compound before the compound is excluded as internally conflicting.
        kind: ``classification`` | ``quantitative`` | ``auto``.
        use_mmp / use_similarity: pairing methods; at least one must be on.
        similarity_threshold: Tanimoto (policy fingerprint) at or above which two compounds are neighbours.
        max_r_atoms: largest varied fragment (heavy atoms) accepted for a matched molecular pair.
        min_outlier_neighbours: a compound is a label outlier when it disagrees with at least this
            many neighbours that all agree with each other.
        keep_pairs: ``all`` | ``cliffs`` (cliffs and undetermined only) — controls ``report.pairs``; counts always cover every pair.
    """
    if not (use_mmp or use_similarity):
        raise ValueError("At least one of use_mmp / use_similarity must be enabled.")
    if keep_pairs not in ("all", "cliffs"):
        raise ValueError("keep_pairs must be 'all' or 'cliffs'")
    pol = policy or ChemicalPolicy.preset("parent")
    vals = [str(r.get(endpoint_field, "")).strip() for r in records if str(r.get(endpoint_field, "")).strip()]
    if kind == "auto":
        kind = "quantitative" if vals and all(CENSOR_RE.match(v) for v in vals) else "classification"
    if kind not in ("classification", "quantitative"):
        raise ValueError("kind must be auto, classification or quantitative")

    # 1. identity groups within context
    groups: Dict[tuple, Dict[str, List[Dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    rep_smiles: Dict[str, str] = {}
    n_records = n_no_value = 0
    for r in records:
        v = str(r.get(endpoint_field, "")).strip()
        smi = r.get("smiles")
        if not smi:
            continue
        if not v:
            n_no_value += 1
            continue
        k = compute_identity(smi)
        if not k.is_valid:
            continue
        key = k.key(pol.identity_level)
        if not key:
            continue
        n_records += 1
        ctx = tuple((c, str(r.get(c, "")).strip()) for c in context_fields)
        groups[ctx][key].append(r)
        rep_smiles.setdefault(key, {"exact": k.exact, "nostereo": k.nostereo, "skeleton": k.nostereo}.get(pol.identity_level, k.parent))

    counters = {"conflict": 0, "unparseable": 0, "missing_units": 0}
    issues: Dict[str, List[Dict[str, Any]]] = {}
    all_pairs: List[CliffPair] = []
    outliers: List[LabelOutlier] = []
    by_rel: Dict[str, int] = defaultdict(int)
    n_comp = n_cliff = n_cons = n_und = 0
    cliff_compounds: set = set()
    trans: Dict[str, Dict[str, Any]] = defaultdict(lambda: {"n": 0, "cliff": 0, "cons": 0, "und": 0, "signed": []})
    scale_seen: Optional[str] = None
    pair_no = 0
    comp_no = 0
    for ctx, by_key in groups.items():
        comps: List[CompoundValue] = []
        for key, members in by_key.items():
            label, lo, hi, disp, n_ex, n_ce, scale = _compound_values(members, endpoint_field, units_field, relation_field, kind, threshold, counters)
            if scale and scale_seen is None:
                scale_seen = scale
            if kind == "quantitative" and lo is None and hi is None:
                continue
            if kind == "classification" and label is None:
                continue
            comp_no += 1
            comps.append(CompoundValue(compound_id=f"m{comp_no:06d}", key=key, smiles=rep_smiles[key], record_ids=[m["record_id"] for m in members],
                                       source_ids=[str(m["source_id"]) for m in members if m.get("source_id") not in (None, "")],
                                       label=label, lo=lo, hi=hi, display=disp, n_exact=n_ex, n_censored=n_ce))
        n_comp += len(comps)
        if len(comps) < 2:
            continue
        found = _find_pairs(comps, pol, use_mmp, use_similarity, similarity_threshold, max_r_atoms)
        context = dict(ctx)
        neighbours: Dict[int, List[Tuple[int, str, str]]] = defaultdict(list)   # i -> [(j, verdict, relationship)]
        for (i, j), d in sorted(found.items()):
            a, b = comps[i], comps[j]
            verdict, mn, mx, signed = _verdict(a, b, kind, threshold)
            pair_no += 1
            by_rel[d["relationship"]] += 1
            if verdict == "cliff":
                n_cliff += 1
                cliff_compounds.update((a.compound_id, b.compound_id))
            elif verdict == "consistent":
                n_cons += 1
            else:
                n_und += 1
            neighbours[i].append((j, verdict, d["relationship"]))
            neighbours[j].append((i, verdict, d["relationship"]))
            if d["transformation"]:
                ra, rb = d["transformation"].split(">>")
                if ra <= rb:
                    tkey, tsigned = f"{ra}>>{rb}", signed
                else:
                    tkey, tsigned = f"{rb}>>{ra}", (-signed if signed is not None else None)
                t = trans[tkey]
                t["n"] += 1
                t[{"cliff": "cliff", "consistent": "cons", "undetermined": "und"}[verdict]] += 1
                if tsigned is not None:
                    t["signed"].append(tsigned)
            if verdict == "cliff":
                for c in (a, b):
                    other = b if c is a else a
                    for rid in c.record_ids:
                        issues.setdefault(rid, []).append({"code": "ACTIVITY_CLIFF", "pair_id": f"p{pair_no:06d}", "relationship": d["relationship"],
                                                           "neighbour_record_ids": other.record_ids, **({"transformation": d["transformation"]} if d["transformation"] else {}),
                                                           **({"min_difference": round(mn, 4)} if mn is not None else {})})
            if keep_pairs == "all" or verdict != "consistent":
                all_pairs.append(CliffPair(pair_id=f"p{pair_no:06d}", context=context, compound_a=a.compound_id, compound_b=b.compound_id,
                                           record_ids_a=a.record_ids, record_ids_b=b.record_ids, source_ids_a=a.source_ids, source_ids_b=b.source_ids, smiles_a=a.smiles, smiles_b=b.smiles,
                                           relationship=d["relationship"], similarity=d["similarity"], core=d["core"], transformation=d["transformation"],
                                           n_changed_atoms=d["n_changed_atoms"], value_a=a.display, value_b=b.display,
                                           min_difference=round(mn, 4) if mn is not None else None, max_difference=round(mx, 4) if mx is not None else None,
                                           signed_difference=round(signed, 4) if signed is not None else None, verdict=verdict))
        # 2. label outliers: every neighbour disagrees with c; no two neighbours disagree with each
        #    other; and at least min_outlier_neighbours of them are proven consistent with one another
        #    (censored neighbours may be undetermined among themselves without blocking the call).
        for i, nb in neighbours.items():
            if len(nb) < min_outlier_neighbours or any(v != "cliff" for _, v, _ in nb):
                continue
            js = [j for j, _, _ in nb]
            verdicts = {}
            blocked = False
            for x in range(len(js)):
                for y in range(x + 1, len(js)):
                    v = _verdict(comps[js[x]], comps[js[y]], kind, threshold)[0]
                    verdicts[(js[x], js[y])] = verdicts[(js[y], js[x])] = v
                    if v == "cliff":
                        blocked = True
            if blocked:
                continue
            n_agree = _largest_consistent_subset(js, verdicts)
            if n_agree < min_outlier_neighbours:
                continue
            c = comps[i]
            outliers.append(LabelOutlier(compound_id=c.compound_id, record_ids=c.record_ids, source_ids=c.source_ids, smiles=c.smiles, value=c.display, context=context,
                                         n_neighbours=len(js), n_agreeing=n_agree, neighbour_ids=[comps[j].compound_id for j in js], neighbour_values=[comps[j].display for j in js],
                                         neighbour_relationships=[rel for _, _, rel in nb]))
            for rid in c.record_ids:
                issues.setdefault(rid, []).append({"code": "LABEL_OUTLIER", "n_neighbours": len(js), "neighbour_record_ids": [r for j in js for r in comps[j].record_ids]})

    tsum = [TransformationSummary(transformation=k, n_pairs=v["n"], n_cliffs=v["cliff"], n_consistent=v["cons"], n_undetermined=v["und"],
                                  mean_signed_difference=round(statistics.fmean(v["signed"]), 4) if v["signed"] else None, n_exact_pairs=len(v["signed"]))
            for k, v in trans.items()]
    tsum.sort(key=lambda t: (-t.n_cliffs, -t.n_pairs, t.transformation))
    decided = n_cliff + n_cons
    return CliffReport(endpoint_field=endpoint_field, context_fields=list(context_fields), identity_level=pol.identity_level, kind=kind, threshold=threshold,
                       scale=scale_seen if kind == "quantitative" else None, use_mmp=use_mmp, use_similarity=use_similarity,
                       similarity_threshold=similarity_threshold, fingerprint=pol.fingerprint, max_r_atoms=max_r_atoms,
                       n_records=n_records, n_no_value=n_no_value, n_compounds=n_comp, n_internal_conflict_excluded=counters["conflict"], n_unparseable=counters["unparseable"],
                       n_missing_units=counters["missing_units"], n_pairs=pair_no, n_pairs_by_relationship=dict(by_rel),
                       n_cliffs=n_cliff, n_consistent=n_cons, n_undetermined=n_und, cliff_fraction=round(n_cliff / decided, 4) if decided else None,
                       n_compounds_in_cliffs=len(cliff_compounds), n_outliers=len(outliers), pairs=all_pairs, transformations=tsum, outliers=outliers,
                       record_issues=issues)


def cliffs_from_annotations(anns: Iterable, policy: ChemicalPolicy, roles: Dict[str, str], threshold: float = 1.0, **kw) -> CliffReport:
    """Adapter for dataset annotations (as produced by the dataset audit)."""
    ctx = [roles[c] for c in ("target",) if c in roles]
    recs = []
    for a in anns:
        if not a.identity or not a.identity.is_valid:
            continue
        d = dict(a.fields)
        d.update({"record_id": a.record_id, "source_id": a.source_id, "smiles": a.parsed_smiles})
        recs.append(d)
    return activity_cliffs(recs, policy, endpoint_field=roles["endpoint"], units_field=roles.get("units"), relation_field=roles.get("relation"),
                           context_fields=ctx, threshold=threshold, **kw)


__all__ = ["activity_cliffs", "cliffs_from_annotations", "mmp_fragments", "CliffReport", "CliffPair", "CompoundValue", "TransformationSummary", "LabelOutlier", "VERDICTS"]
