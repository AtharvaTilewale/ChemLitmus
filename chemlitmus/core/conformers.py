"""Multiple conformer generation via RDKit ETKDG."""
from __future__ import annotations
from typing import Optional, List
from pydantic import BaseModel

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

class ConformerResult(BaseModel):
    input_smiles: str
    num_generated: int = 0
    error: Optional[str] = None

def generate_conformers(smiles: str, num_conformers: int = 50, output_sdf: Optional[str] = None) -> ConformerResult:
    """Generate multiple 3D conformers for a given SMILES."""
    if not _RDKIT_AVAILABLE:
        return ConformerResult(input_smiles=smiles, error="RDKit is not installed.")
        
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ConformerResult(input_smiles=smiles, error="Invalid SMILES string.")

    try:
        # Add hydrogens for 3D generation
        mol = Chem.AddHs(mol)
        
        # Embed using ETKDG default parameters
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, randomSeed=42)
        
        # Optimize with MMFF
        try:
            AllChem.MMFFOptimizeMoleculeConfs(mol)
        except Exception:
            pass # fallback to unoptimized if MMFF fails

        generated_count = len(conf_ids)
        
        if output_sdf and generated_count > 0:
            writer = Chem.SDWriter(output_sdf)
            for cid in conf_ids:
                writer.write(mol, confId=cid)
            writer.close()
            
        return ConformerResult(input_smiles=smiles, num_generated=generated_count)
    except Exception as e:
        return ConformerResult(input_smiles=smiles, error=str(e))
