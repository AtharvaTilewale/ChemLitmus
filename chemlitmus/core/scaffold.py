"""Murcko Scaffold extraction."""
from __future__ import annotations
from typing import Optional
from pydantic import BaseModel

try:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

class ScaffoldResult(BaseModel):
    input_smiles: str
    scaffold_smiles: Optional[str] = None
    error: Optional[str] = None
    
    @property
    def success(self) -> bool:
        return self.error is None

def extract_scaffold(smiles: str) -> ScaffoldResult:
    """Extract the Murcko framework from a SMILES string."""
    if not _RDKIT_AVAILABLE:
        return ScaffoldResult(input_smiles=smiles, error="RDKit is not installed.")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ScaffoldResult(input_smiles=smiles, error="Invalid SMILES string.")

    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smi = Chem.MolToSmiles(core)
        return ScaffoldResult(input_smiles=smiles, scaffold_smiles=scaffold_smi)
    except Exception as e:
        return ScaffoldResult(input_smiles=smiles, error=str(e))
