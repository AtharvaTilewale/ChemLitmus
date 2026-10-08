"""ChemLitmus: High-performance SMILES validation and PubChem lookup."""

from pathlib import Path
from typing import Optional, Union, List, Callable, Dict
import concurrent.futures

from chemlitmus.config import settings
from chemlitmus.core.database import DatabaseManager
from chemlitmus.core.pubchem import PubChemClient, PubChemCompound
from chemlitmus.core.smiles import SMILESValidationResult, validate_smiles, generate_structure
from chemlitmus.core.standardize import (
    StandardizeResult,
    StepResult,
    standardize_smiles,
    VALID_STEPS as STANDARDIZE_STEPS,
)
from chemlitmus.core.tautomers import (
    TautomerResult,
    enumerate_tautomers,
)
from chemlitmus.core.rgroup import rgroup_decomposition, RGroupResult
from chemlitmus.core.augment import augment_smiles, AugmentResult
from chemlitmus.core.atommap import map_atoms, AtomMapResult
from chemlitmus.core.reaction import validate_reaction, ReactionResult
from chemlitmus.core.conformers import generate_conformers, ConformerResult
from chemlitmus.core.scaffold import extract_scaffold, ScaffoldResult
from chemlitmus.core.stereo import analyze_stereochemistry, StereoResult
from chemlitmus.core.iupacname import (
    IUPACResult,
    get_iupac_name,
)
from chemlitmus.core.cheminfo import (
    FingerprintResult,
    FilterResult,
    SimilarityResult,
    SubstructureHit,
    compute_fingerprint,
    apply_filters,
    compute_similarity,
    substructure_search,
)
from chemlitmus.core.identity import (
    IDENTITY_LEVELS,
    IdentityKeys,
    IdentityGroup,
    IdentityReport,
    compute_identity,
    group_by_identity,
    strictest_shared_level,
    describe_difference,
)
from chemlitmus.core.libdiff import DiffEntry, LibraryDiff, diff_libraries
from chemlitmus.core.diagnose import (
    DIAGNOSTIC_CATEGORIES,
    SmilesProblem,
    SmilesDiagnosis,
    diagnose_smiles,
)
from chemlitmus.core.smartsaudit import (
    SmartsAuditResult,
    PatternAudit,
    SensitivitySummary,
    SmartsExplanation,
    AtomExplanation,
    PREPARATIONS,
    AUDIT_CHECKS,
    audit_smarts,
    explain_smarts,
    load_patterns,
    load_reference_library,
    prepare_molecule,
)
from chemlitmus.utils.parsers import parse_compounds_file
from chemlitmus.utils.export import export_results

__version__ = "1.0.0"
__author__ = "Atharva Tilewale"
__license__ = "MIT"


def lookup(
    query: Union[str, int],
    search_type: str = "auto",
    use_cache: bool = True,
) -> Optional[PubChemCompound]:
    """Lookup chemical compound information from PubChem with local caching."""
    query_str = str(query).strip()
    db = DatabaseManager()

    if use_cache:
        db.init_db()
        cached = db.get_cached_compound(query_str)
        if cached:
            return cached

    client = PubChemClient()
    compound = client.lookup(query_str, search_type=search_type)

    if compound and use_cache:
        db.cache_compound(query_str, compound)
        if compound.canonical_smiles:
            db.cache_compound(compound.canonical_smiles, compound)
        if compound.cid:
            db.cache_compound(str(compound.cid), compound)

    return compound

def lookup_by_name(name: str, use_cache: bool = True) -> Optional[PubChemCompound]:
    """
    Search PubChem for a chemical by its common or IUPAC name.
    
    Args:
        name: The chemical name (e.g., "Aspirin", "benzene")
        use_cache: Whether to use the local SQLite database
        
    Returns:
        PubChemCompound object or None if not found
    """
    name_str = str(name).strip()
    db = DatabaseManager() if use_cache else None
    
    if use_cache:
        db.init_db()
        cached = db.get_cached_compound(name_str)
        if cached:
            return cached
            
    # Fetch from PubChem
    client = PubChemClient()
    result = client.lookup_by_name(name_str)
    
    if result and use_cache and result:
        db.cache_compound(name_str, result)
        if result.canonical_smiles:
            db.cache_compound(result.canonical_smiles, result)
        if result.cid:
            db.cache_compound(str(result.cid), result)
        
    return result

def lookup_file(
    input_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    output_format: str = "csv",
    remove_duplicates: bool = True,
    progress_callback: Optional[Callable] = None,
) -> List[PubChemCompound]:
    """Process a batch file of SMILES/CIDs using multithreading."""
    input_path = Path(input_file)
    queries = parse_compounds_file(input_path)
    
    if remove_duplicates:
        queries = list(dict.fromkeys(queries))
        
    results = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=settings.max_workers) as executor:
        # Submit all tasks to the thread pool
        future_to_query = {executor.submit(lookup, query): query for query in queries}
        
        for future in concurrent.futures.as_completed(future_to_query):
            query = future_to_query[future]
            try:
                compound = future.result()
                if compound:
                    results.append(compound)
                else:
                    results.append(PubChemCompound(input_query=query))
            except Exception:
                results.append(PubChemCompound(input_query=query))
            
            if progress_callback:
                progress_callback()

    if output_file:
        export_results(results, Path(output_file), output_format)
        
    return results


def download_structure(
    cid: int,
    format: str = "sdf",
    dimension: str = "3d",
    output_dir: str = "structures",
    force: bool = False,
) -> str:
    """Download a single chemical structure from PubChem."""
    client = PubChemClient()
    return client.download_structure(cid, format, dimension, output_dir, force)

__all__ = [
    # Molecular identity and library comparison
    "IDENTITY_LEVELS",
    "IdentityKeys",
    "IdentityGroup",
    "IdentityReport",
    "compute_identity",
    "group_by_identity",
    "strictest_shared_level",
    "describe_difference",
    "DiffEntry",
    "LibraryDiff",
    "diff_libraries",
    # SMILES diagnostics
    "DIAGNOSTIC_CATEGORIES",
    "SmilesProblem",
    "SmilesDiagnosis",
    "diagnose_smiles",
    # SMARTS auditing
    "SmartsAuditResult",
    "PatternAudit",
    "SensitivitySummary",
    "SmartsExplanation",
    "AtomExplanation",
    "PREPARATIONS",
    "AUDIT_CHECKS",
    "audit_smarts",
    "explain_smarts",
    "load_patterns",
    "load_reference_library",
    "prepare_molecule",
    "lookup",
    "lookup_by_name",
    "lookup_file",
    "download_structure",
    "generate_structure",
    "validate_smiles",
    "compute_fingerprint",
    "apply_filters",
    "compute_similarity",
    "substructure_search",
    "SubstructureHit",
    "standardize_smiles",
    "StandardizeResult",
    "STANDARDIZE_STEPS",
    "get_iupac_name",
    "IUPACResult",
    "enumerate_tautomers",
    "TautomerResult",
    "rgroup_decomposition",
    "RGroupResult",
    "augment_smiles",
    "AugmentResult",
    "map_atoms",
    "AtomMapResult",
    "validate_reaction",
    "ReactionResult",
    "generate_conformers",
    "ConformerResult",
    "extract_scaffold",
    "ScaffoldResult",
    "analyze_stereochemistry",
    "StereoResult",
    "SMILESValidationResult",
    "FingerprintResult",
    "FilterResult",
    "SimilarityResult",
    "PubChemClient",
    "PubChemCompound",
    "DatabaseManager",
]
