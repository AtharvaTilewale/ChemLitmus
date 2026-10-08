"""Reaction SMILES validation and analysis."""
from __future__ import annotations
from typing import Optional, List
from pydantic import BaseModel

try:
    from rdkit.Chem import rdChemReactions
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

class ReactionResult(BaseModel):
    """Result model for reaction analysis."""
    input_smiles: str
    is_valid: bool = False
    num_reactants: int = 0
    num_products: int = 0
    num_agents: int = 0
    error: Optional[str] = None

def validate_reaction(smiles: str) -> ReactionResult:
    """Validate a Reaction SMILES (SMIRKS) string."""
    if not _RDKIT_AVAILABLE:
        return ReactionResult(input_smiles=smiles, error="RDKit is not installed.")

    try:
        rxn = rdChemReactions.ReactionFromSmarts(smiles, useSmiles=True)
        if rxn is None:
            return ReactionResult(input_smiles=smiles, error="Could not parse Reaction SMILES.")
        
        return ReactionResult(
            input_smiles=smiles,
            is_valid=True,
            num_reactants=rxn.GetNumReactantTemplates(),
            num_products=rxn.GetNumProductTemplates(),
            num_agents=rxn.GetNumAgentTemplates(),
        )
    except Exception as e:
        return ReactionResult(input_smiles=smiles, error=str(e))
