"""Cheminformatics utilities: fingerprints, similarity search, and drug-likeness filters.

All functions are fully offline and powered by RDKit. No network access is required.
"""

from typing import List, Optional, Union

from pydantic import BaseModel, Field

try:
    from rdkit import Chem
    from rdkit.Chem import DataStructs, Descriptors, FilterCatalog, MACCSkeys, QED, rdMolDescriptors
    from rdkit.Chem import rdFingerprintGenerator

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


class FingerprintResult(BaseModel):
    """Result model for a single molecular fingerprint computation."""

    smiles: str
    fingerprint_type: str
    n_bits: int
    n_on_bits: int
    density: float
    bit_string: str
    hex_string: str


class RuleResult(BaseModel):
    """Pass/fail result for a single drug-likeness rule."""

    passed: bool
    details: str = ""


class FilterResult(BaseModel):
    """Result model for drug-likeness and ADMET filter evaluation."""

    smiles: str
    molecular_weight: Optional[float] = None
    logp: Optional[float] = None
    hbd: Optional[int] = None
    hba: Optional[int] = None
    tpsa: Optional[float] = None
    rotatable_bonds: Optional[int] = None
    heavy_atom_count: Optional[int] = None
    molar_refractivity: Optional[float] = None
    qed_score: Optional[float] = None
    lipinski: Optional[RuleResult] = None
    veber: Optional[RuleResult] = None
    ghose: Optional[RuleResult] = None
    egan: Optional[RuleResult] = None
    ro3: Optional[RuleResult] = None
    pains: Optional[RuleResult] = None
    passes_all: bool = False
    error: Optional[str] = None


class SimilarityResult(BaseModel):
    """Result model for a single similarity search hit."""

    rank: int
    query: str
    hit: str
    similarity: float
    fingerprint_type: str


_VALID_FP_TYPES = {"ecfp4", "ecfp6", "fcfp4", "maccs", "rdkit", "atompair", "torsion"}




class SubstructureHit(BaseModel):
    """Result model for a single substructure match."""

    smiles: str
    matched: bool
    match_indices: List[int] = Field(default_factory=list)

def _get_fp_generator(fp_type: str, n_bits: int = 2048):
    """Return an RDKit fingerprint generator for the given type."""
    t = fp_type.lower()
    if t == "ecfp4":
        return rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n_bits)
    elif t == "ecfp6":
        return rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=n_bits)
    elif t == "fcfp4":
        return rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n_bits)
    elif t == "rdkit":
        return rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=n_bits)
    elif t == "atompair":
        return rdFingerprintGenerator.GetAtomPairGenerator(fpSize=n_bits)
    elif t == "torsion":
        return rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=n_bits)
    elif t == "maccs":
        return None
    else:
        raise ValueError(f"Unknown fingerprint type: {fp_type!r}")


def _compute_single_fp(mol, fp_type: str, n_bits: int = 2048):
    """Compute a fingerprint for an RDKit mol. Returns (bit_vector, actual_bits)."""
    t = fp_type.lower()
    if t == "maccs":
        fp = MACCSkeys.GenMACCSKeys(mol)
        return fp, 167
    gen = _get_fp_generator(t, n_bits)
    fp = gen.GetFingerprint(mol)
    return fp, n_bits


def _fp_to_hex(fp) -> str:
    """Convert an RDKit bit vector to a compact hex string."""
    bit_str = fp.ToBitString()
    pad = (4 - len(bit_str) % 4) % 4
    padded = "0" * pad + bit_str
    return "".join(format(int(padded[i : i + 4], 2), "x") for i in range(0, len(padded), 4))


def compute_fingerprint(
    smiles: str,
    fp_type: str = "ecfp4",
    n_bits: int = 2048,
) -> Union[FingerprintResult, List[FingerprintResult]]:
    """Compute molecular fingerprint(s) from a SMILES string (offline, RDKit-based).

    Args:
        smiles:  Input SMILES string.
        fp_type: Fingerprint algorithm. One of: ecfp4, ecfp6, fcfp4, maccs, rdkit,
                 atompair, torsion, or 'all' to compute all types at once.
        n_bits:  Number of bits for hashed fingerprints (ignored for MACCS, fixed at 167).

    Returns:
        A single FingerprintResult, or a list of FingerprintResult if fp_type='all'.
    """
    if not _RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required. Install with: pip install rdkit")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles!r}")

    def _build(t: str) -> FingerprintResult:
        fp, actual_bits = _compute_single_fp(mol, t, n_bits)
        on = fp.GetNumOnBits()
        return FingerprintResult(
            smiles=smiles,
            fingerprint_type=t,
            n_bits=actual_bits,
            n_on_bits=on,
            density=round(on / actual_bits, 6) if actual_bits else 0.0,
            bit_string=fp.ToBitString(),
            hex_string=_fp_to_hex(fp),
        )

    if fp_type.lower() == "all":
        return [_build(t) for t in sorted(_VALID_FP_TYPES)]
    if fp_type.lower() not in _VALID_FP_TYPES:
        raise ValueError(
            f"Unknown fp_type: {fp_type!r}. Valid: {sorted(_VALID_FP_TYPES | {'all'})}"
        )
    return _build(fp_type.lower())


def apply_filters(
    smiles: str,
    rules: Optional[List[str]] = None,
) -> FilterResult:
    """Evaluate drug-likeness and ADMET filters on a SMILES string (offline, RDKit-based).

    Args:
        smiles: Input SMILES string.
        rules:  Rule names to apply: lipinski, veber, ghose, egan, ro3, pains, qed, or all.
                Pass None or ['all'] to apply every rule.

    Returns:
        A FilterResult with per-rule pass/fail and computed property values.
    """
    if not _RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required. Install with: pip install rdkit")

    _valid_rules = {"lipinski", "veber", "ghose", "egan", "ro3", "pains", "qed"}
    if rules is None or (len(rules) == 1 and rules[0].lower() == "all"):
        active_rules = _valid_rules
    else:
        active_rules = set()
        for r in rules:
            rc = r.strip().lower()
            if rc not in _valid_rules:
                raise ValueError(
                    f"Unknown rule: {rc!r}. Valid options: {sorted(_valid_rules | {'all'})}"
                )
            active_rules.add(rc)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return FilterResult(smiles=smiles, error="Invalid SMILES: could not be parsed by RDKit.")

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = Descriptors.TPSA(mol)
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    heavy = mol.GetNumHeavyAtoms()
    mr = Descriptors.MolMR(mol)
    qed_val = QED.qed(mol)

    result = FilterResult(
        smiles=smiles,
        molecular_weight=round(mw, 3),
        logp=round(logp, 3),
        hbd=hbd,
        hba=hba,
        tpsa=round(tpsa, 3),
        rotatable_bonds=rot,
        heavy_atom_count=heavy,
        molar_refractivity=round(mr, 3),
        qed_score=round(qed_val, 4),
    )

    rule_pass_list: List[bool] = []

    if "lipinski" in active_rules:
        v = []
        if mw > 500: v.append(f"MW {mw:.1f}>500")
        if logp > 5: v.append(f"LogP {logp:.2f}>5")
        if hbd > 5: v.append(f"HBD {hbd}>5")
        if hba > 10: v.append(f"HBA {hba}>10")
        p = len(v) == 0
        result.lipinski = RuleResult(
            passed=p,
            details="Passes Lipinski Ro5" if p else f"Violations: {'; '.join(v)}",
        )
        rule_pass_list.append(p)

    if "veber" in active_rules:
        v = []
        if rot > 10: v.append(f"RotBonds {rot}>10")
        if tpsa > 140: v.append(f"TPSA {tpsa:.1f}>140")
        p = len(v) == 0
        result.veber = RuleResult(
            passed=p,
            details="Passes Veber rules" if p else f"Violations: {'; '.join(v)}",
        )
        rule_pass_list.append(p)

    if "ghose" in active_rules:
        v = []
        if not (160 <= mw <= 480): v.append(f"MW {mw:.1f} not in [160,480]")
        if not (-0.4 <= logp <= 5.6): v.append(f"LogP {logp:.2f} not in [-0.4,5.6]")
        if not (20 <= heavy <= 70): v.append(f"HeavyAtoms {heavy} not in [20,70]")
        if not (40 <= mr <= 130): v.append(f"MR {mr:.1f} not in [40,130]")
        p = len(v) == 0
        result.ghose = RuleResult(
            passed=p,
            details="Passes Ghose filter" if p else f"Violations: {'; '.join(v)}",
        )
        rule_pass_list.append(p)

    if "egan" in active_rules:
        v = []
        if tpsa > 131.6: v.append(f"TPSA {tpsa:.1f}>131.6")
        if logp > 5.88: v.append(f"LogP {logp:.2f}>5.88")
        p = len(v) == 0
        result.egan = RuleResult(
            passed=p,
            details="Passes Egan filter" if p else f"Violations: {'; '.join(v)}",
        )
        rule_pass_list.append(p)

    if "ro3" in active_rules:
        v = []
        if mw > 300: v.append(f"MW {mw:.1f}>300")
        if logp > 3: v.append(f"LogP {logp:.2f}>3")
        if hbd > 3: v.append(f"HBD {hbd}>3")
        if hba > 3: v.append(f"HBA {hba}>3")
        p = len(v) == 0
        result.ro3 = RuleResult(
            passed=p,
            details="Passes Ro3 (lead-like)" if p else f"Violations: {'; '.join(v)}",
        )
        rule_pass_list.append(p)

    if "pains" in active_rules:
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)
        entry = catalog.GetFirstMatch(mol)
        has_alert = entry is not None
        result.pains = RuleResult(
            passed=not has_alert,
            details="No PAINS alert" if not has_alert else f"PAINS alert: {entry.GetDescription()}",
        )
        rule_pass_list.append(not has_alert)

    # QED is informational only - not added to the pass list
    result.passes_all = all(rule_pass_list) if rule_pass_list else False
    return result


def compute_similarity(
    query_smiles: str,
    library: List[str],
    fp_type: str = "ecfp4",
    n_bits: int = 2048,
    threshold: float = 0.0,
    top_n: Optional[int] = None,
) -> List[SimilarityResult]:
    """Compute Tanimoto similarity between a query SMILES and a list of library SMILES (offline).

    Args:
        query_smiles: Query SMILES string.
        library:      List of library SMILES strings.
        fp_type:      Fingerprint type (ecfp4, ecfp6, fcfp4, maccs, rdkit, atompair, torsion).
        n_bits:       Fingerprint size in bits.
        threshold:    Minimum Tanimoto similarity (0.0-1.0) to include in results.
        top_n:        Maximum number of results. None returns all above threshold.

    Returns:
        List of SimilarityResult sorted by similarity descending.
    """
    if not _RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required. Install with: pip install rdkit")

    fp_type = fp_type.lower()
    if fp_type not in _VALID_FP_TYPES:
        raise ValueError(f"Unknown fp_type: {fp_type!r}. Valid: {sorted(_VALID_FP_TYPES)}")

    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        raise ValueError(f"Invalid query SMILES: {query_smiles!r}")

    query_fp, _ = _compute_single_fp(query_mol, fp_type, n_bits)
    results: List[SimilarityResult] = []

    for lib_smi in library:
        if not lib_smi or not lib_smi.strip():
            continue
        lib_mol = Chem.MolFromSmiles(lib_smi.strip())
        if lib_mol is None:
            continue
        lib_fp, _ = _compute_single_fp(lib_mol, fp_type, n_bits)
        sim = DataStructs.TanimotoSimilarity(query_fp, lib_fp)
        if sim >= threshold:
            results.append(
                SimilarityResult(
                    rank=0,
                    query=query_smiles,
                    hit=lib_smi.strip(),
                    similarity=round(sim, 6),
                    fingerprint_type=fp_type,
                )
            )

    results.sort(key=lambda x: x.similarity, reverse=True)
    if top_n is not None and top_n > 0:
        results = results[:top_n]
    for i, r in enumerate(results, 1):
        r.rank = i
    return results


def substructure_search(
    query: str,
    library: List[str],
    is_smarts: bool = True,
) -> List[SubstructureHit]:
    """Search a library for compounds containing a specific substructure.

    Args:
        query: SMILES or SMARTS string representing the substructure.
        library: List of SMILES strings to search against.
        is_smarts: If True, parse the query as SMARTS. If False, parse as SMILES.

    Returns:
        A list of SubstructureHit objects for compounds that matched the query.
    """
    if not _RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is not installed.")

    # Parse the query
    if is_smarts:
        q_mol = Chem.MolFromSmarts(query)
        if q_mol is None:
            raise ValueError(f"Invalid SMARTS query: {query}")
    else:
        q_mol = Chem.MolFromSmiles(query)
        if q_mol is None:
            raise ValueError(f"Invalid SMILES query: {query}")

    hits = []
    for smi in library:
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        if mol.HasSubstructMatch(q_mol):
            # Get the indices of the atoms in the target molecule that match the query
            match = mol.GetSubstructMatch(q_mol)
            hits.append(
                SubstructureHit(
                    smiles=smi,
                    matched=True,
                    match_indices=list(match),
                )
            )

    return hits

