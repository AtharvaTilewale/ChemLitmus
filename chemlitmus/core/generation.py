"""Evaluate a collection of generated molecules.

Denominators are explicit and reported with every metric:

* **validity** = parseable / all generated attempts (empty strings and duplicates included);
* **uniqueness** = distinct keys / valid outputs, at each identity level separately;
* **novelty** = outputs with no match in the supplied reference set(s) / *unique valid* outputs,
  and never claimed beyond the references actually checked.

Repaired candidates are evaluated separately from the raw outputs: a repair makes a string parse,
it does not make it the molecule the model intended. Descriptor compliance and alert prevalence
are reported as properties, never as evidence of potency or synthesisability.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Dict, List, Optional, Sequence

from pydantic import BaseModel, Field

from chemlitmus.core.identity import IDENTITY_LEVELS, compute_identity
from chemlitmus.core.policy import ChemicalPolicy

try:
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.DisableLog("rdApp.*")
except ImportError:  # pragma: no cover
    Chem = None


class MoleculeEval(BaseModel):
    index: int
    input: str
    valid: bool
    smiles: Optional[str] = None
    repaired: bool = False
    repair_candidate: Optional[str] = None
    duplicate_of: Optional[int] = Field(None, description="Index of the first output with the same key at the policy identity level.")
    novel: Optional[bool] = Field(None, description="None when no reference set was supplied.")
    nearest_reference: Optional[str] = None
    nearest_similarity: Optional[float] = None
    scaffold: Optional[str] = None
    alerts: List[str] = Field(default_factory=list)
    constraints_failed: List[str] = Field(default_factory=list)
    descriptors: Dict[str, float] = Field(default_factory=dict)
    reason: Optional[str] = Field(None, description="Why an output is invalid.")


class GenerationReport(BaseModel):
    n_generated: int = Field(description="All attempts, including empty strings.")
    n_valid: int
    validity: float = Field(description="n_valid / n_generated.")
    n_empty: int = 0
    invalid_reasons: Dict[str, int] = Field(default_factory=dict)
    identity_level: str
    n_unique: int = Field(description="Distinct keys at the policy identity level among valid outputs.")
    uniqueness: float = Field(description="n_unique / n_valid (undefined -> 0.0 with n_valid = 0; see 'undefined' list).")
    uniqueness_by_level: Dict[str, float] = Field(default_factory=dict, description="Per identity level; levels nest, so these are not independent.")
    reference_sets: List[str] = Field(default_factory=list, description="Names of the reference collections actually checked.")
    n_reference: int = 0
    n_novel: Optional[int] = None
    novelty: Optional[float] = Field(None, description="n_novel / n_unique (declared population: unique valid outputs).")
    novelty_by_level: Dict[str, float] = Field(default_factory=dict)
    n_scaffolds: int = 0
    scaffold_diversity: Optional[float] = Field(None, description="distinct scaffolds / unique valid outputs.")
    nearest_neighbour_similarity: Dict[str, float] = Field(default_factory=dict, description="min / median / max over valid outputs with a reference set.")
    fingerprint: Optional[str] = None
    descriptor_summary: Dict[str, Dict[str, float]] = Field(default_factory=dict)
    alerts_by_set: Dict[str, int] = Field(default_factory=dict)
    alert_denominator: int = 0
    constraints: Dict[str, int] = Field(default_factory=dict, description="constraint -> number of valid outputs satisfying it.")
    n_repaired: int = 0
    repaired_report: Optional["GenerationReport"] = Field(None, description="The same metrics computed over repaired candidates only.")
    undefined: List[str] = Field(default_factory=list, description="Metrics that are undefined for this input (e.g. no valid outputs), rather than reported as 0.")
    molecules: List[MoleculeEval] = Field(default_factory=list)
    note: str = ("Validity is over all attempts; uniqueness over valid outputs; novelty over unique valid outputs against the "
                 "reference sets listed. Descriptor compliance and alert counts describe structure, not potency, safety or synthesisability.")


def _fp(mol, policy: ChemicalPolicy):
    from chemlitmus.core.cheminfo import _compute_single_fp
    return _compute_single_fp(mol, policy.fingerprint, policy.fingerprint_bits)[0]


def _descriptors(mol) -> Dict[str, float]:
    from chemlitmus.core.dataset_audit import _descriptors as d
    return d(mol)


def _check_constraints(desc: Dict[str, float], constraints: Dict[str, Any]) -> List[str]:
    failed = []
    for name, spec in constraints.items():
        v = desc.get(name)
        if v is None:
            failed.append(f"{name}: not computed")
            continue
        lo, hi = (spec if isinstance(spec, (list, tuple)) else (None, spec))
        if lo is not None and v < lo:
            failed.append(f"{name}={v} < {lo}")
        if hi is not None and v > hi:
            failed.append(f"{name}={v} > {hi}")
    return failed


def evaluate_generated(
    generated: Sequence[str],
    reference: Optional[Dict[str, Sequence[str]]] = None,
    policy: Optional[ChemicalPolicy] = None,
    constraints: Optional[Dict[str, Any]] = None,
    repair: bool = False,
    keep_molecules: bool = True,
) -> GenerationReport:
    """Evaluate generated SMILES.

    Args:
        generated: raw model outputs, exactly as produced (do not pre-filter).
        reference: ``{name: [smiles, ...]}`` collections to test novelty against.
        constraints: ``{descriptor: max}`` or ``{descriptor: (min, max)}`` using the descriptor
            names of the audit (mw, logp, hbd, hba, tpsa, rotatable_bonds, heavy_atoms, rings, fraction_csp3).
        repair: also evaluate mechanical repair candidates for invalid outputs, separately.
    """
    pol = policy or ChemicalPolicy.preset("parent")
    mols: List[MoleculeEval] = []
    reasons = Counter()
    n_empty = 0
    for i, s in enumerate(generated):
        txt = "" if s is None else str(s)
        me = MoleculeEval(index=i, input=txt, valid=False)
        if not txt.strip():
            n_empty += 1
            me.reason = "empty"
            reasons["empty"] += 1
            mols.append(me)
            continue
        from chemlitmus.core.smiles import SmilesParseError, mol_from_smiles
        try:
            mol = mol_from_smiles(txt)
        except SmilesParseError:
            mol, me.reason = None, "whitespace"
        if mol is None:
            me.reason = me.reason or "unparseable"
            reasons[me.reason] += 1
            if repair:
                from chemlitmus.core.diagnose import diagnose_smiles
                d = diagnose_smiles(txt, try_repair=True)
                if d.repair_status == "candidate":
                    me.repair_candidate = d.repaired_smiles
            mols.append(me)
            continue
        me.valid = True
        me.smiles = Chem.MolToSmiles(mol)
        mols.append(me)

    valid = [m for m in mols if m.valid]
    keys = {m.index: compute_identity(m.smiles) for m in valid}
    seen: Dict[str, int] = {}
    for m in valid:
        k = keys[m.index].key(pol.identity_level)
        if k in seen:
            m.duplicate_of = seen[k]
        else:
            seen[k] = m.index
    unique = [m for m in valid if m.duplicate_of is None]

    undefined: List[str] = []
    if not valid:
        undefined += ["uniqueness", "novelty", "scaffold_diversity", "nearest_neighbour_similarity"]
    elif not unique:  # pragma: no cover - cannot happen, kept for clarity
        undefined.append("novelty")

    rep = GenerationReport(
        n_generated=len(generated), n_valid=len(valid), validity=round(len(valid) / len(generated), 4) if generated else 0.0,
        n_empty=n_empty, invalid_reasons=dict(reasons), identity_level=pol.identity_level,
        n_unique=len(unique), uniqueness=round(len(unique) / len(valid), 4) if valid else 0.0, undefined=undefined,
    )
    rep.uniqueness_by_level = {lvl: round(len({keys[m.index].key(lvl) for m in valid if keys[m.index].key(lvl)}) / len(valid), 4) for lvl in IDENTITY_LEVELS} if valid else {}

    # descriptors, scaffolds, alerts, constraints over valid outputs
    from chemlitmus.core.dataset_audit import _alerts
    for m in valid:
        mol = Chem.MolFromSmiles(m.smiles)
        m.descriptors = _descriptors(mol)
        try:
            sc = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
        except Exception:
            sc = ""
        m.scaffold = sc or "acyclic"
        m.alerts = [a["description"] for a in _alerts(mol, pol) if a["set"]]
        if constraints:
            m.constraints_failed = _check_constraints(m.descriptors, constraints)
    if valid:
        for name in ["mw", "logp", "tpsa", "heavy_atoms", "rings", "fraction_csp3"]:
            vals = sorted(m.descriptors[name] for m in valid if name in m.descriptors)
            if vals:
                rep.descriptor_summary[name] = {"n": float(len(vals)), "min": vals[0], "median": vals[len(vals) // 2], "max": vals[-1],
                                                "mean": round(sum(vals) / len(vals), 4)}
        alert_counter = Counter()
        for m in valid:
            if m.alerts:
                for a in {x for x in m.alerts}:
                    alert_counter[a] += 1
        rep.alerts_by_set = {k: v for k, v in alert_counter.most_common(20)}
        rep.alert_denominator = len(valid)
        scafs = {m.scaffold for m in unique if m.scaffold != "acyclic"}
        rep.n_scaffolds = len(scafs)
        rep.scaffold_diversity = round(len({m.scaffold for m in unique}) / len(unique), 4) if unique else None
        if constraints:
            rep.constraints = {name: sum(1 for m in valid if not any(f.startswith(name) for f in m.constraints_failed)) for name in constraints}

    # novelty against references
    if reference:
        ref_keys: Dict[str, set] = {lvl: set() for lvl in IDENTITY_LEVELS}
        ref_fps, ref_smiles = [], []
        for name, smis in reference.items():
            for s in smis:
                k = compute_identity(s)
                if not k.is_valid:
                    continue
                for lvl in IDENTITY_LEVELS:
                    if k.key(lvl):
                        ref_keys[lvl].add(k.key(lvl))
                mol = Chem.MolFromSmiles(s)
                if mol is not None:
                    ref_fps.append(_fp(mol, pol)); ref_smiles.append(s)
        rep.reference_sets = list(reference)
        rep.n_reference = len(ref_smiles)
        for m in valid:
            k = keys[m.index].key(pol.identity_level)
            m.novel = k not in ref_keys[pol.identity_level]
        if unique:
            rep.n_novel = sum(1 for m in unique if m.novel)
            rep.novelty = round(rep.n_novel / len(unique), 4)
            rep.novelty_by_level = {lvl: round(sum(1 for m in unique if keys[m.index].key(lvl) not in ref_keys[lvl]) / len(unique), 4) for lvl in IDENTITY_LEVELS}
        if ref_fps and valid:
            sims = []
            for m in valid:
                mol = Chem.MolFromSmiles(m.smiles)
                s = DataStructs.BulkTanimotoSimilarity(_fp(mol, pol), ref_fps)
                j = max(range(len(s)), key=s.__getitem__)
                m.nearest_similarity, m.nearest_reference = round(s[j], 4), ref_smiles[j]
                sims.append(m.nearest_similarity)
            sims.sort()
            rep.nearest_neighbour_similarity = {"min": sims[0], "median": sims[len(sims) // 2], "max": sims[-1]}
            rep.fingerprint = f"{pol.fingerprint}/{pol.fingerprint_bits}"

    # repaired candidates, evaluated as their own population
    cands = [m.repair_candidate for m in mols if m.repair_candidate]
    rep.n_repaired = len(cands)
    if repair and cands:
        sub = evaluate_generated(cands, reference=reference, policy=pol, constraints=constraints, repair=False, keep_molecules=False)
        for m in sub.molecules:
            m.repaired = True
        rep.repaired_report = sub
    if keep_molecules:
        rep.molecules = mols
    return rep


GenerationReport.model_rebuild()

__all__ = ["evaluate_generated", "GenerationReport", "MoleculeEval"]
