"""Core module exposing SMILES validation, PubChem client, and database manager."""

from chemlitmus.core.database import DatabaseManager
from chemlitmus.core.pubchem import PubChemClient, PubChemCompound
from chemlitmus.core.smiles import SMILESValidationResult, validate_smiles
from chemlitmus.core.identity import (
    IdentityKeys,
    IdentityGroup,
    IdentityReport,
    compute_identity,
    group_by_identity,
)
from chemlitmus.core.libdiff import DiffEntry, LibraryDiff, diff_libraries
from chemlitmus.core.diagnose import SmilesProblem, SmilesDiagnosis, diagnose_smiles
from chemlitmus.core.smartsaudit import (
    SmartsAuditResult,
    PatternAudit,
    SensitivitySummary,
    SmartsExplanation,
    audit_smarts,
    explain_smarts,
    load_patterns,
    load_reference_library,
)

__all__ = [
    "validate_smiles",
    "SMILESValidationResult",
    "PubChemClient",
    "PubChemCompound",
    "DatabaseManager",
    "IdentityKeys",
    "IdentityGroup",
    "IdentityReport",
    "compute_identity",
    "group_by_identity",
    "DiffEntry",
    "LibraryDiff",
    "diff_libraries",
    "SmilesProblem",
    "SmilesDiagnosis",
    "diagnose_smiles",
    "SmartsAuditResult",
    "PatternAudit",
    "SensitivitySummary",
    "SmartsExplanation",
    "audit_smarts",
    "explain_smarts",
    "load_patterns",
    "load_reference_library",
]