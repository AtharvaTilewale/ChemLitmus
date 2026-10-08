"""Stereochemistry analysis."""
from __future__ import annotations
from typing import Optional, List, Dict, Any
from pydantic import BaseModel

try:
    from rdkit import Chem
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

class StereoResult(BaseModel):
    input_smiles: str
    chiral_centers: List[Dict[str, Any]] = [] # e.g. [{"atom_idx": 1, "config": "R"}]
    has_unassigned: bool = False
    error: Optional[str] = None
    
    @property
    def success(self) -> bool:
        return self.error is None

def analyze_stereochemistry(smiles: str) -> StereoResult:
    """Analyze stereocenters and unassigned centers in a SMILES."""
    if not _RDKIT_AVAILABLE:
        return StereoResult(input_smiles=smiles, error="RDKit is not installed.")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return StereoResult(input_smiles=smiles, error="Invalid SMILES string.")

    try:
        Chem.AssignStereochemistry(mol, force=True, flagPossibleStereoCenters=True)
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        
        chiral_list = []
        unassigned = False
        for idx, config in centers:
            chiral_list.append({"atom_idx": idx, "config": config})
            if config == "?":
                unassigned = True
                
        return StereoResult(
            input_smiles=smiles, 
            chiral_centers=chiral_list,
            has_unassigned=unassigned,
        )
    except Exception as e:
        return StereoResult(input_smiles=smiles, error=str(e))
