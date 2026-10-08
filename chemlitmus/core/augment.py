"""SMILES Augmentation (Randomization)."""
from __future__ import annotations
from typing import Optional, List
from pydantic import BaseModel
import random

try:
    from rdkit import Chem
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

class AugmentResult(BaseModel):
    input_smiles: str
    augmented_smiles: List[str] = []
    error: Optional[str] = None
    
    @property
    def success(self) -> bool:
        return self.error is None

def augment_smiles(smiles: str, num_augmentations: int = 5) -> AugmentResult:
    """Generate multiple uncanonical randomized SMILES strings for data augmentation."""
    if not _RDKIT_AVAILABLE:
        return AugmentResult(input_smiles=smiles, error="RDKit is not installed.")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return AugmentResult(input_smiles=smiles, error="Invalid SMILES string.")

    try:
        # Keep track of unique SMILES
        unique_smiles = set()
        
        # We might need many attempts to get N unique SMILES if the molecule is small
        max_attempts = num_augmentations * 10
        attempts = 0
        
        while len(unique_smiles) < num_augmentations and attempts < max_attempts:
            # doRandom=True generates randomized SMILES
            random_smi = Chem.MolToSmiles(mol, doRandom=True, canonical=False)
            unique_smiles.add(random_smi)
            attempts += 1
            
        return AugmentResult(
            input_smiles=smiles,
            augmented_smiles=list(unique_smiles)
        )
    except Exception as e:
        return AugmentResult(input_smiles=smiles, error=str(e))
