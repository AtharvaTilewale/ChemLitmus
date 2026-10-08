"""Core module exposing SMILES validation, PubChem client, and database manager."""

from chemlitmus.core.database import DatabaseManager
from chemlitmus.core.pubchem import PubChemClient, PubChemCompound
from chemlitmus.core.smiles import SMILESValidationResult, validate_smiles
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
    "SmartsAuditResult",
    "PatternAudit",
    "SensitivitySummary",
    "SmartsExplanation",
    "audit_smarts",
    "explain_smarts",
    "load_patterns",
    "load_reference_library",
]