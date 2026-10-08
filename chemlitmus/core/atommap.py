"""Atom Mapping."""
from __future__ import annotations
from typing import Optional
from pydantic import BaseModel

try:
    from rdkit import Chem
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

class AtomMapResult(BaseModel):
    input_smiles: str
    mapped_smiles: Optional[str] = None
    error: Optional[str] = None
    
    @property
    def success(self) -> bool:
        return self.error is None

def map_atoms(smiles: str) -> AtomMapResult:
    """Assign unique atom map numbers to all atoms in a molecule."""
    if not _RDKIT_AVAILABLE:
        return AtomMapResult(input_smiles=smiles, error="RDKit is not installed.")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return AtomMapResult(input_smiles=smiles, error="Invalid SMILES string.")

    try:
        # Clear existing mappings just in case
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
            
        # Assign new mappings (1-indexed)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
            
        mapped_smi = Chem.MolToSmiles(mol)
        return AtomMapResult(input_smiles=smiles, mapped_smiles=mapped_smi)
    except Exception as e:
        return AtomMapResult(input_smiles=smiles, error=str(e))
