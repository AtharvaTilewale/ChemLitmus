"""Unified dataset audit.

One pass over a record set that answers: what is wrong or uncertain in this chemical dataset,
how does it affect my workflow, and what is the evidence? Every finding is an :class:`Issue`
with a stable code, a severity, the affected records and evidence; every record keeps its
original fields; every count reconciles with the input.

The audit runs offline. It never contacts an external service.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from statistics import median
from typing import Any, Dict, List, Optional, Sequence, Tuple

from pydantic import BaseModel, Field

from chemlitmus.core.identity import IDENTITY_LEVELS, IdentityKeys, _differs_by, compute_identity, describe_difference
from chemlitmus.core.policy import ChemicalPolicy, resolve_policy
from chemlitmus.core.records import Issue, Record, RecordSet, read_records

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, FilterCatalog, rdMolDescriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.DisableLog("rdApp.*")
except ImportError:  # pragma: no cover
    Chem = None

AUDIT_SCHEMA_VERSION = "1"

# ------------------------------------------------------------------------------------ issue catalogue

ISSUE_CATALOGUE: Dict[str, Tuple[str, str, str]] = {
    # code: (default severity, description, suggested action)
    "PARSE_EMPTY": ("warning", "Structure field is empty.", "Supply a structure or drop the record explicitly."),
    "PARSE_WHITESPACE": ("error", "Structure field contains whitespace; RDKit would silently keep only the first token.", "Fix the field; if it is a SMILES-plus-name line, split it."),
    "PARSE_INVALID": ("error", "RDKit cannot parse the structure.", "Review the candidate repair (if any) or re-export from the source."),
    "SDF_PARSE_FAILED": ("error", "RDKit could not read or sanitise the SDF record.", "Inspect the MOL block; check valences and bond blocks."),
    "REPAIR_CANDIDATE": ("info", "A mechanical repair produces a parseable structure.", "Verify against the source before accepting; apply only under repair.mode=apply."),
    "STD_CHANGED_TEXT": ("info", "Standardisation changed the SMILES text but not the molecular identity.", "No action; recorded for provenance."),
    "STD_CHANGED_IDENTITY": ("warning", "Standardisation changed the molecular identity (fragment, charge, tautomer or stereo).", "Confirm the policy intends this; the original is retained."),
    "STD_FAILED": ("error", "A standardisation step raised.", "Record kept with the parsed structure; inspect manually."),
    "FRAG_MULTIPLE": ("info", "The record has more than one component.", "Check that largest-fragment selection picked the intended component."),
    "FRAG_MULTIPLE_ORGANIC": ("warning", "More than one substantial organic component (mixture, co-crystal or unresolved salt pair).", "Review: the parent chosen by size may not be the active compound."),
    "CHARGE_NET": ("info", "The parent carries a net charge after neutralisation (permanent charge).", "No action unless the dataset expects neutral species."),
    "ISOTOPE": ("info", "Isotopic labels present.", "Decide whether isotopes are meaningful for this dataset (policy.isotopes)."),
    "STEREO_UNASSIGNED": ("info", "Stereocentres present without assignment.", "Treat as racemic/unknown; identity at 'nostereo' will merge with assigned forms."),
    "STEREO_PARTIAL": ("warning", "Some stereocentres assigned, others not.", "Check whether the record is a partially defined mixture."),
    "ELEMENT_UNUSUAL": ("info", "Elements outside the common organic set.", "Organometallics and inorganics may be outside the scope of descriptors and alerts."),
    "DESC_MW_HIGH": ("info", "Molecular weight above 1000 Da.", "Confirm the record is a small molecule."),
    "DESC_TINY": ("info", "Fewer than 3 heavy atoms.", "Likely a solvent, ion or fragment; check the source."),
    "DUP_EXACT": ("warning", "Exact duplicate of an earlier record.", "Deduplicate or confirm independent measurements."),
    "DUP_PARENT": ("warning", "Same parent compound as an earlier record (salt, counter-ion or charge form).", "Group or deduplicate according to the identity policy."),
    "DUP_TAUTOMER": ("info", "Tautomer of an earlier record.", "Group if tautomers are one compound for this purpose."),
    "DUP_NOSTEREO": ("info", "Stereoisomer (or stereo-unassigned form) of an earlier record.", "Keep distinct unless stereo is irrelevant."),
    "DUP_SKELETON": ("info", "Same 2D skeleton as an earlier record (differs in stereo and/or tautomer).", "Review."),
    "ALERT_MATCH": ("info", "Matches a structural alert.", "Context, not a verdict: alerts flag assay-interference motifs, not toxicity or activity."),
    "LABEL_CONFLICT": ("warning", "Records of the same compound carry contradictory labels or measurements.", "Review; do not resolve by keeping the strongest value."),
    "ACTIVITY_CLIFF": ("info", "A near-identical compound (matched molecular pair or high-similarity neighbour) carries a proven label difference at or above the threshold.", "SAR when both values are sound, an error when one is not; the pair is evidence for review, not a verdict."),
    "LABEL_OUTLIER": ("warning", "Every near neighbour disagrees with this compound's value while agreeing with each other.", "Check the source record: this is the pattern a wrong value produces most often, but genuine cliffs produce it too."),
    "UNITS_MISSING": ("warning", "A quantitative measurement has no units.", "Supply units; the value cannot be compared."),
    "RELATION_CENSORED": ("info", "Measurement is censored (<, >, range).", "Treat as a bound, not a value."),
    "SPLIT_OVERLAP": ("error", "A compound appears in more than one split (at the policy identity level).", "Move or drop one occurrence; re-run the leakage report."),
    "SPLIT_RELATED": ("info", "A test-set compound is closely related to a training compound (shared scaffold or high similarity).", "Relatedness is an evaluation-design consideration, not a defect."),
    "SPLIT_MISSING": ("warning", "Record has no split label.", "Assign a split or exclude from leakage analysis."),
    "TEMPORAL_ORDER": ("warning", "A test record is dated no later than the latest training record.", "Review the temporal split definition."),
}
"""Every issue code the audit can emit: (default severity, meaning, suggested action)."""

SEVERITY_RANK = {"ignore": 0, "info": 1, "warning": 2, "error": 3}
COMMON_ELEMENTS = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53}
DESCRIPTOR_NAMES = ["mw", "logp", "hbd", "hba", "tpsa", "rotatable_bonds", "heavy_atoms", "rings", "fraction_csp3"]


def make_issue(code: str, policy: ChemicalPolicy, message: Optional[str] = None, field: Optional[str] = None, **evidence) -> Optional[Issue]:
    sev, desc, _ = ISSUE_CATALOGUE[code]
    sev = policy.severity.overrides.get(code, sev)
    if sev == "ignore":
        return None
    return Issue(code=code, severity=sev, message=message or desc, field=field, evidence=evidence)


# ------------------------------------------------------------------------------------ models

class Transformation(BaseModel):
    """One standardisation step with before/after and whether identity changed."""

    operation: str
    parameters: Dict[str, Any] = Field(default_factory=dict)
    before: str
    after: str
    changed_text: bool
    changed_identity: bool = Field(description="Canonical SMILES differ (exact-level identity changed).")
    identity_relation: Optional[str] = Field(None, description="How before and after relate (identity-level language) when identity changed.")
    error: Optional[str] = None


class RecordAnnotation(BaseModel):
    record_id: str
    position: int
    source_id: Optional[str] = None
    status: str
    input_structure: Optional[str] = None
    parsed_smiles: Optional[str] = None
    standardized_smiles: Optional[str] = None
    transformations: List[Transformation] = Field(default_factory=list)
    standardization_changed_identity: bool = False
    repair_candidate: Optional[str] = None
    repair_status: Optional[str] = None
    identity: Optional[IdentityKeys] = None
    group_id: Optional[str] = Field(None, description="Identity group at the policy level.")
    n_components: int = 0
    n_organic_components: int = 0
    net_charge: Optional[int] = None
    has_isotopes: bool = False
    n_stereocentres: int = 0
    n_unassigned_stereocentres: int = 0
    descriptors: Dict[str, float] = Field(default_factory=dict)
    scaffold: Optional[str] = Field(None, description="Murcko scaffold SMILES; 'acyclic' when none.")
    alerts: List[Dict[str, Any]] = Field(default_factory=list, description="[{set, description, atoms}]")
    split: Optional[str] = None
    issues: List[Issue] = Field(default_factory=list)
    fields: Dict[str, Any] = Field(default_factory=dict)

    @property
    def max_severity(self) -> str:
        ranks = [SEVERITY_RANK[i.severity] for i in self.issues] or [0]
        return {v: k for k, v in SEVERITY_RANK.items()}[max(ranks)]


class IdentityGroupOut(BaseModel):
    group_id: str
    level: str
    key: str
    size: int
    record_ids: List[str]
    source_ids: List[Optional[str]]
    differs_by: List[str] = Field(default_factory=list)
    splits: List[str] = Field(default_factory=list, description="Distinct split labels among members (leakage when > 1).")


class IssueRow(BaseModel):
    code: str
    severity: str
    record_id: Optional[str]
    source_id: Optional[str] = None
    message: str
    evidence: Dict[str, Any] = Field(default_factory=dict)
    suggested_action: str = ""


class DescriptorSummary(BaseModel):
    n: int
    min: float
    median: float
    max: float
    mean: float


class AuditSummary(BaseModel):
    n_total: int
    n_ok: int
    n_empty: int
    n_invalid: int
    n_unsupported: int
    n_error: int
    n_annotated: int = Field(description="Records with a parsed structure that went through standardisation and identity.")
    processing_complete: bool = Field(description="False when any computation failed for a record (STD_FAILED) or a preparation/match failed.")
    issues_by_code: Dict[str, int] = Field(default_factory=dict)
    issues_by_severity: Dict[str, int] = Field(default_factory=dict)
    records_by_max_severity: Dict[str, int] = Field(default_factory=dict)
    identity_level: str
    n_groups: int = 0
    n_collapsed: int = Field(0, description="Records sharing a key with an earlier record at the policy level.")
    n_groups_by_level: Dict[str, int] = Field(default_factory=dict)
    descriptors: Dict[str, DescriptorSummary] = Field(default_factory=dict)
    alerts_by_set: Dict[str, int] = Field(default_factory=dict, description="Records matching >= 1 alert of each set (denominator: n_annotated).")
    alert_denominator: int = 0
    dataset_warnings: List[str] = Field(default_factory=list, description="Configuration-level findings (not per-record).")


class DatasetAudit(BaseModel):
    schema_version: str = AUDIT_SCHEMA_VERSION
    source_file: str
    policy: ChemicalPolicy
    policy_hash: str
    roles: Dict[str, str]
    summary: AuditSummary
    records: List[RecordAnnotation]
    groups: List[IdentityGroupOut]
    issues: List[IssueRow]
    leakage: Optional[Dict[str, Any]] = Field(None, description="LeakageReport (when a split column is present).")
    label_conflicts: Optional[Dict[str, Any]] = Field(None, description="LabelConflictReport (when an endpoint column is present).")


# ------------------------------------------------------------------------------------ per-record chemistry

def _standardize(mol, policy: ChemicalPolicy) -> Tuple[Optional["Chem.Mol"], List[Transformation]]:
    steps: List[Transformation] = []
    cur = mol

    def _step(op: str, fn, **params):
        nonlocal cur
        before = Chem.MolToSmiles(cur)
        try:
            new = fn(cur)
            after = Chem.MolToSmiles(new)
        except Exception as exc:
            steps.append(Transformation(operation=op, parameters=params, before=before, after=before, changed_text=False, changed_identity=False, error=f"{type(exc).__name__}: {exc}"))
            return
        rel = None
        if after != before:
            ka, kb = compute_identity(before), compute_identity(after)
            rel = describe_difference(ka, kb) if ka.is_valid and kb.is_valid else "unknown"
        steps.append(Transformation(operation=op, parameters=params, before=before, after=after, changed_text=after != before,
                                    changed_identity=after != before, identity_relation=rel))
        cur = new

    if policy.fragment == "largest_organic":
        _step("largest_organic_fragment", rdMolStandardize.LargestFragmentChooser().choose, chooser="rdkit LargestFragmentChooser")
    if policy.neutralize:
        _step("neutralize", rdMolStandardize.Uncharger().uncharge, method="rdkit Uncharger")
    if policy.isotopes == "strip":
        def _strip_iso(m):
            m2 = Chem.Mol(m)
            for a in m2.GetAtoms():
                a.SetIsotope(0)
            return m2
        _step("strip_isotopes", _strip_iso)
    if policy.stereo == "strip":
        def _strip_stereo(m):
            m2 = Chem.Mol(m); Chem.RemoveStereochemistry(m2); return m2
        _step("strip_stereo", _strip_stereo)
    if policy.tautomer == "canonicalize":
        _step("canonical_tautomer", rdMolStandardize.TautomerEnumerator().Canonicalize, method="rdkit TautomerEnumerator")
    return cur, steps


def _components(mol, policy: ChemicalPolicy) -> Tuple[int, int]:
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    organic = sum(1 for f in frags if f.GetNumHeavyAtoms() >= policy.organic_component_min_heavy_atoms and any(a.GetAtomicNum() == 6 for a in f.GetAtoms()))
    return len(frags), organic


def _stereo_counts(mol) -> Tuple[int, int]:
    centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    return len(centres), sum(1 for _, tag in centres if tag == "?")


def _descriptors(mol) -> Dict[str, float]:
    return {
        "mw": round(Descriptors.MolWt(mol), 3), "logp": round(Descriptors.MolLogP(mol), 3),
        "hbd": float(rdMolDescriptors.CalcNumHBD(mol)), "hba": float(rdMolDescriptors.CalcNumHBA(mol)),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 3), "rotatable_bonds": float(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "heavy_atoms": float(mol.GetNumHeavyAtoms()), "rings": float(rdMolDescriptors.CalcNumRings(mol)),
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
    }


def _scaffold(mol) -> str:
    try:
        sc = MurckoScaffold.GetScaffoldForMol(mol)
        s = Chem.MolToSmiles(sc) if sc is not None else ""
        return s or "acyclic"
    except Exception:
        return "acyclic"


_CATALOGS: Dict[Tuple[str, ...], Any] = {}


def _catalog(sets: Sequence[str]):
    key = tuple(sets)
    if key not in _CATALOGS:
        params = FilterCatalog.FilterCatalogParams()
        for s in sets:
            params.AddCatalog(getattr(FilterCatalog.FilterCatalogParams.FilterCatalogs, s))
        _CATALOGS[key] = FilterCatalog.FilterCatalog(params)
    return _CATALOGS[key]


def _alerts(mol, policy: ChemicalPolicy) -> List[Dict[str, Any]]:
    if not policy.alert_sets:
        return []
    from chemlitmus.core.smartsaudit import prepare_molecule_status
    pm, status = prepare_molecule_status(mol, policy.preparation)
    if pm is None:
        return [{"set": "", "description": f"not evaluated: {status}", "atoms": []}]
    out = []
    for entry in _catalog(policy.alert_sets).GetMatches(pm):
        atoms: List[int] = []
        for fm in entry.GetFilterMatches(pm):
            atoms.extend(sorted({j for _, j in fm.atomPairs}))
        out.append({"set": entry.GetProp("FilterSet") if "FilterSet" in list(entry.GetPropList()) else "", "description": entry.GetDescription(), "atoms": sorted(set(atoms))})
    return out


def annotate_record(rec: Record, policy: ChemicalPolicy, split: Optional[str] = None) -> RecordAnnotation:
    ann = RecordAnnotation(record_id=rec.record_id, position=rec.position, source_id=rec.source_id, status=rec.status,
                           input_structure=rec.structure if rec.structure_format == "smiles" else "(molblock)",
                           parsed_smiles=rec.parsed_smiles, issues=list(rec.issues), fields=dict(rec.fields), split=split)
    for iss in ann.issues:  # apply severity overrides to ingestion issues
        iss.severity = policy.severity.overrides.get(iss.code, iss.severity)
    ann.issues = [i for i in ann.issues if i.severity != "ignore"]
    if rec.status != "ok" or not rec.parsed_smiles:
        if rec.status == "invalid" and policy.repair.mode != "none" and rec.structure_format == "smiles" and rec.structure:
            from chemlitmus.core.diagnose import diagnose_smiles
            d = diagnose_smiles(str(rec.structure), try_repair=True)
            ann.repair_status = d.repair_status
            if d.repair_status == "candidate":
                ann.repair_candidate = d.repaired_smiles
                iss = make_issue("REPAIR_CANDIDATE", policy, candidate=d.repaired_smiles, edits=d.repairs_applied)
                if iss:
                    ann.issues.append(iss)
        return ann
    mol = Chem.MolFromSmiles(rec.parsed_smiles)
    # components / charge / isotopes / stereo / elements on the *input* molecule
    ann.n_components, ann.n_organic_components = _components(mol, policy)
    if ann.n_components > 1:
        iss = make_issue("FRAG_MULTIPLE", policy, n_components=ann.n_components)
        if iss:
            ann.issues.append(iss)
    if policy.flag_multiple_organic_components and ann.n_organic_components > 1:
        iss = make_issue("FRAG_MULTIPLE_ORGANIC", policy, n_organic_components=ann.n_organic_components)
        if iss:
            ann.issues.append(iss)
    ann.has_isotopes = any(a.GetIsotope() for a in mol.GetAtoms())
    if ann.has_isotopes:
        iss = make_issue("ISOTOPE", policy)
        if iss:
            ann.issues.append(iss)
    unusual = sorted({a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() not in COMMON_ELEMENTS})
    if unusual:
        iss = make_issue("ELEMENT_UNUSUAL", policy, elements=unusual)
        if iss:
            ann.issues.append(iss)
    # standardisation
    std, steps = _standardize(mol, policy)
    ann.transformations = steps
    for s in steps:
        if s.error:
            iss = make_issue("STD_FAILED", policy, operation=s.operation, error=s.error)
            if iss:
                ann.issues.append(iss)
    ann.standardized_smiles = Chem.MolToSmiles(std)
    if ann.standardized_smiles != rec.parsed_smiles:
        rel = next((s.identity_relation for s in steps if s.changed_identity), None)
        ann.standardization_changed_identity = True
        iss = make_issue("STD_CHANGED_IDENTITY", policy, before=rec.parsed_smiles, after=ann.standardized_smiles, relation=rel)
        if iss:
            ann.issues.append(iss)
    ann.net_charge = Chem.GetFormalCharge(std)
    if ann.net_charge:
        iss = make_issue("CHARGE_NET", policy, charge=ann.net_charge)
        if iss:
            ann.issues.append(iss)
    ann.n_stereocentres, ann.n_unassigned_stereocentres = _stereo_counts(std)
    if ann.n_unassigned_stereocentres:
        code = "STEREO_PARTIAL" if ann.n_unassigned_stereocentres < ann.n_stereocentres else "STEREO_UNASSIGNED"
        iss = make_issue(code, policy, n_centres=ann.n_stereocentres, n_unassigned=ann.n_unassigned_stereocentres)
        if iss:
            ann.issues.append(iss)
    ann.descriptors = _descriptors(std)
    if ann.descriptors["mw"] > 1000:
        iss = make_issue("DESC_MW_HIGH", policy, mw=ann.descriptors["mw"])
        if iss:
            ann.issues.append(iss)
    if ann.descriptors["heavy_atoms"] < 3:
        iss = make_issue("DESC_TINY", policy, heavy_atoms=ann.descriptors["heavy_atoms"])
        if iss:
            ann.issues.append(iss)
    ann.scaffold = _scaffold(std)
    ann.identity = compute_identity(rec.parsed_smiles)
    for a in _alerts(std, policy):
        ann.alerts.append(a)
        if a["set"] or a["description"].startswith("not evaluated"):
            iss = make_issue("ALERT_MATCH", policy, message=f"{a['set']}: {a['description']}" if a["set"] else a["description"], **a)
            if iss:
                ann.issues.append(iss)
    return ann


# ------------------------------------------------------------------------------------ dataset level

def _group(anns: List[RecordAnnotation], level: str, policy: ChemicalPolicy) -> List[IdentityGroupOut]:
    by_key: Dict[str, List[RecordAnnotation]] = defaultdict(list)
    for a in anns:
        if a.identity and a.identity.is_valid and a.identity.key(level):
            by_key[a.identity.key(level)].append(a)
    groups: List[IdentityGroupOut] = []
    for n, (key, members) in enumerate(sorted(by_key.items(), key=lambda kv: (-len(kv[1]), kv[1][0].position))):
        gid = f"g{n:06d}"
        for m in members:
            m.group_id = gid
        groups.append(IdentityGroupOut(group_id=gid, level=level, key=key, size=len(members), record_ids=[m.record_id for m in members],
                                       source_ids=[m.source_id for m in members],
                                       differs_by=_differs_by([m.identity for m in members], level) if len(members) > 1 else [],
                                       splits=sorted({m.split for m in members if m.split})))
    # duplicate issues: every member after the first, at each level separately (not additive)
    for lvl, code in (("exact", "DUP_EXACT"), ("parent", "DUP_PARENT"), ("tautomer", "DUP_TAUTOMER"), ("nostereo", "DUP_NOSTEREO"), ("skeleton", "DUP_SKELETON")):
        seen: Dict[str, RecordAnnotation] = {}
        for a in sorted(anns, key=lambda x: x.position):
            k = a.identity.key(lvl) if a.identity and a.identity.is_valid else None
            if not k:
                continue
            if k in seen:
                first = seen[k]
                # only report the strictest applicable level per pair: skip looser codes when a stricter one already fired for this pair
                stricter = [i for i in a.issues if i.code.startswith("DUP_") and i.evidence.get("first_record_id") == first.record_id]
                if stricter:
                    continue
                iss = make_issue(code, policy, level=lvl, first_record_id=first.record_id, first_source_id=first.source_id,
                                 relation=describe_difference(first.identity, a.identity))
                if iss:
                    a.issues.append(iss)
            else:
                seen[k] = a
    return groups


def _descriptor_summary(anns: List[RecordAnnotation]) -> Dict[str, DescriptorSummary]:
    out = {}
    for name in DESCRIPTOR_NAMES:
        vals = [a.descriptors[name] for a in anns if name in a.descriptors]
        if vals:
            out[name] = DescriptorSummary(n=len(vals), min=min(vals), median=float(median(vals)), max=max(vals), mean=round(sum(vals) / len(vals), 4))
    return out


def audit_dataset(
    source: str | RecordSet,
    policy: Optional[str | dict | ChemicalPolicy] = None,
    structure_column: Optional[str] = None,
    id_column: Optional[str] = None,
    roles: Optional[Dict[str, str]] = None,
    ambiguous: str = "error",
    progress_callback=None,
) -> DatasetAudit:
    """Audit a dataset file (or an already-read :class:`RecordSet`) under a chemical policy."""
    pol = resolve_policy(policy)
    rs = source if isinstance(source, RecordSet) else read_records(source, structure_column=structure_column, id_column=id_column, roles=roles, ambiguous=ambiguous)
    split_col = rs.roles.get("split")
    anns: List[RecordAnnotation] = []
    for n, rec in enumerate(rs.records, 1):
        split = (str(rec.fields.get(split_col)).strip() or None) if split_col and rec.fields.get(split_col) not in (None, "") else None
        anns.append(annotate_record(rec, pol, split=split))
        if progress_callback:
            progress_callback(n, len(rs.records))
    if split_col:
        for a in anns:
            if a.status == "ok" and a.split is None:
                iss = make_issue("SPLIT_MISSING", pol, field=split_col)
                if iss:
                    a.issues.append(iss)
    groups = _group([a for a in anns if a.identity], pol.identity_level, pol)
    n_groups_by_level = {lvl: len({a.identity.key(lvl) for a in anns if a.identity and a.identity.is_valid and a.identity.key(lvl)}) for lvl in IDENTITY_LEVELS}
    annotated = [a for a in anns if a.identity and a.identity.is_valid]

    leak = labels = None
    if split_col and annotated:
        from chemlitmus.core.leakage import leakage_from_annotations
        leak = leakage_from_annotations(annotated, pol, groups, date_field=rs.roles.get("date"))
        for a in annotated:
            for i in leak.record_issues.get(a.record_id, []):
                iss = make_issue(i["code"], pol, **{k: v for k, v in i.items() if k != "code"})
                if iss:
                    a.issues.append(iss)
    if rs.roles.get("endpoint") and annotated:
        from chemlitmus.core.labels import conflicts_from_annotations
        labels = conflicts_from_annotations(annotated, pol, rs.roles)
        for a in annotated:
            for i in labels.record_issues.get(a.record_id, []):
                iss = make_issue(i["code"], pol, **{k: v for k, v in i.items() if k != "code"})
                if iss:
                    a.issues.append(iss)

    issue_rows: List[IssueRow] = []
    for a in anns:
        for i in a.issues:
            issue_rows.append(IssueRow(code=i.code, severity=i.severity, record_id=a.record_id, source_id=a.source_id, message=i.message,
                                       evidence=i.evidence, suggested_action=ISSUE_CATALOGUE.get(i.code, ("", "", ""))[2]))
    by_sev = Counter(i.severity for i in issue_rows)
    alert_sets = Counter()
    for a in annotated:
        for s in {al["set"] for al in a.alerts if al["set"]}:
            alert_sets[s] += 1
    warnings: List[str] = list(rs.notes)
    if split_col is None:
        warnings.append("no split column: leakage checks not run")
    if rs.roles.get("endpoint") is None:
        warnings.append("no endpoint column: label-conflict checks not run")
    processing_complete = not any(i.code == "STD_FAILED" for i in issue_rows) and not any(al["description"].startswith("not evaluated") for a in annotated for al in a.alerts)
    summary = AuditSummary(
        n_total=rs.n_total, n_ok=rs.n_ok, n_empty=rs.n_empty, n_invalid=rs.n_invalid, n_unsupported=rs.n_unsupported, n_error=rs.n_error,
        n_annotated=len(annotated), processing_complete=processing_complete,
        issues_by_code=dict(Counter(i.code for i in issue_rows)), issues_by_severity=dict(by_sev),
        records_by_max_severity=dict(Counter(a.max_severity for a in anns)),
        identity_level=pol.identity_level, n_groups=len(groups), n_collapsed=sum(g.size - 1 for g in groups), n_groups_by_level=n_groups_by_level,
        descriptors=_descriptor_summary(annotated), alerts_by_set=dict(alert_sets), alert_denominator=len(annotated), dataset_warnings=warnings,
    )
    return DatasetAudit(source_file=rs.source_file, policy=pol, policy_hash=pol.hash, roles=rs.roles, summary=summary, records=anns, groups=groups,
                        issues=issue_rows, leakage=leak.model_dump() if leak else None, label_conflicts=labels.model_dump() if labels else None)


__all__ = ["audit_dataset", "annotate_record", "DatasetAudit", "AuditSummary", "RecordAnnotation", "IdentityGroupOut", "IssueRow",
           "Transformation", "ISSUE_CATALOGUE", "AUDIT_SCHEMA_VERSION", "make_issue", "SEVERITY_RANK"]
