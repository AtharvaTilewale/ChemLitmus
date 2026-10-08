"""Tautomer enumeration using RDKit MolStandardize.

Critical for protein-ligand docking preparation as different tautomers 
can bind with vastly different affinities.
"""

from __future__ import annotations

from typing import List, Optional
from pydantic import BaseModel, Field

try:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


class TautomerResult(BaseModel):
    """Result model for tautomer enumeration."""

    input_smiles: str
    tautomers: List[str] = Field(default_factory=list)
    num_tautomers: int = 0
    canonical_tautomer: Optional[str] = None
    error: Optional[str] = None

    @property
    def success(self) -> bool:
        return self.error is None


def enumerate_tautomers(
    smiles: str,
    max_tautomers: int = 1000,
) -> TautomerResult:
    """Enumerate all plausible tautomers for a given SMILES string.

    Args:
        smiles: Input SMILES string.
        max_tautomers: Maximum number of tautomers to generate (default: 1000).

    Returns:
        TautomerResult containing the generated tautomer SMILES strings and
        the canonical tautomer form.
    """
    if not _RDKIT_AVAILABLE:
        return TautomerResult(
            input_smiles=smiles,
            error="RDKit is not installed. Install it with: pip install rdkit",
        )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return TautomerResult(
            input_smiles=smiles,
            error=f"Invalid SMILES: could not parse '{smiles}'",
        )

    try:
        enumerator = rdMolStandardize.TautomerEnumerator()
        enumerator.SetMaxTautomers(max_tautomers)
        
        # Enumerate all tautomers
        tauts = enumerator.Enumerate(mol)
        tautomer_smiles = []
        for t in tauts:
            smi = Chem.MolToSmiles(t)
            if smi not in tautomer_smiles:
                tautomer_smiles.append(smi)
                
        # Find the canonical tautomer
        canon_mol = enumerator.Canonicalize(mol)
        canon_smiles = Chem.MolToSmiles(canon_mol)
        
        return TautomerResult(
            input_smiles=smiles,
            tautomers=tautomer_smiles,
            num_tautomers=len(tautomer_smiles),
            canonical_tautomer=canon_smiles,
        )
    except Exception as exc:
        return TautomerResult(
            input_smiles=smiles,
            error=f"Tautomer enumeration failed: {exc}",
        )
