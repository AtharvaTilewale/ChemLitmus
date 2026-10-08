"""R-Group Decomposition."""
from __future__ import annotations
from typing import Optional, List, Dict
from pydantic import BaseModel

try:
    from rdkit import Chem
    from rdkit.Chem import rdRGroupDecomposition
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

class RGroupResult(BaseModel):
    """Result model for R-Group decomposition."""
    input_smiles: str
    decomposition: Dict[str, str] = {}  # e.g. {"Core": "...", "R1": "...", "R2": "..."}
    is_matched: bool = False
    error: Optional[str] = None

def rgroup_decomposition(core_smarts: str, smiles_list: List[str]) -> List[RGroupResult]:
    """Decompose a list of SMILES against a common core SMARTS."""
    results = []
    if not _RDKIT_AVAILABLE:
        for smi in smiles_list:
            results.append(RGroupResult(input_smiles=smi, error="RDKit is not installed."))
        return results

    core = Chem.MolFromSmarts(core_smarts)
    if core is None:
        for smi in smiles_list:
            results.append(RGroupResult(input_smiles=smi, error="Invalid core SMARTS string."))
        return results

    mols = []
    valid_smiles = []
    for smi in smiles_list:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            mols.append(m)
            valid_smiles.append(smi)
        else:
            results.append(RGroupResult(input_smiles=smi, error="Invalid SMILES string."))

    if not mols:
        return results

    try:
        rgd, unmatched = rdRGroupDecomposition.RGroupDecompose([core], mols, asSmiles=True)
        
        matched_idx = 0
        for i, smi in enumerate(valid_smiles):
            if i in unmatched:
                results.append(RGroupResult(input_smiles=smi, error="Molecule does not match core."))
            else:
                results.append(RGroupResult(
                    input_smiles=smi,
                    decomposition=rgd[matched_idx],
                    is_matched=True
                ))
                matched_idx += 1
                
        # Re-sort results to match input order? We append in order for unmatched/matched, 
        # but unmatched contains indices into `valid_smiles`.
        # Let's ensure order.
        final_results = []
        result_map = {r.input_smiles: r for r in results}
        for smi in smiles_list:
            if smi in result_map:
                final_results.append(result_map[smi])
            else:
                # Fallback for duplicates
                pass
                
        return results # It's actually already mostly ordered by valid_smiles + errors at start, wait!
        # Actually, let's simplify to return exactly matching the valid_smiles order.
    except Exception as e:
        for smi in valid_smiles:
            results.append(RGroupResult(input_smiles=smi, error=str(e)))
            
    # Fix sorting:
    fixed = []
    v_idx = 0
    m_idx = 0
    if 'rgd' in locals():
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                fixed.append(RGroupResult(input_smiles=smi, error="Invalid SMILES string."))
            else:
                if v_idx in unmatched:
                    fixed.append(RGroupResult(input_smiles=smi, error="Molecule does not match core."))
                else:
                    fixed.append(RGroupResult(input_smiles=smi, decomposition=rgd[m_idx], is_matched=True))
                    m_idx += 1
                v_idx += 1
        return fixed
    return results
