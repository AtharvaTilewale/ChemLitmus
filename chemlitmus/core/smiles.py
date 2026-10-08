"""SMILES validation, property calculation, and 2D/3D structure generation using RDKit."""

from pathlib import Path
from typing import Optional, Union
from pydantic import BaseModel, Field

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
    from rdkit import RDLogger
    RDKIT_AVAILABLE = True
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    RDKIT_AVAILABLE = False


class SMILESValidationResult(BaseModel):
    """Result model for SMILES validation and property extraction."""

    input_smiles: str
    is_valid: bool
    canonical_smiles: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    logp: Optional[float] = None
    hbd: Optional[int] = Field(None, description="H-Bond Donors")
    hba: Optional[int] = Field(None, description="H-Bond Acceptors")
    tpsa: Optional[float] = Field(None, description="Topological Polar Surface Area")
    heavy_atom_count: Optional[int] = None
    error_message: Optional[str] = None


def validate_smiles(smiles_str: str) -> SMILESValidationResult:
    """
    Validate and canonicalize a SMILES string, returning chemical descriptors.

    Args:
        smiles_str: Input SMILES string to validate.

    Returns:
        SMILESValidationResult object.
    """
    if not smiles_str or not isinstance(smiles_str, str):
        return SMILESValidationResult(
            input_smiles=str(smiles_str),
            is_valid=False,
            error_message="Empty or invalid input type",
        )

    smiles_clean = smiles_str.strip()

    if not RDKIT_AVAILABLE:
        return SMILESValidationResult(
            input_smiles=smiles_clean,
            is_valid=True,
            canonical_smiles=smiles_clean,
            error_message="RDKit is not installed; skipping detailed validation",
        )

    mol = Chem.MolFromSmiles(smiles_clean)
    if mol is None:
        return SMILESValidationResult(
            input_smiles=smiles_clean,
            is_valid=False,
            error_message="RDKit failed to parse SMILES string",
        )

    try:
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        mw = float(Descriptors.ExactMolWt(mol))
        formula = rdMolDescriptors.CalcMolFormula(mol)
        logp = float(Descriptors.MolLogP(mol))
        hbd = int(rdMolDescriptors.CalcNumHBD(mol))
        hba = int(rdMolDescriptors.CalcNumHBA(mol))
        tpsa = float(Descriptors.TPSA(mol))
        heavy_atoms = int(mol.GetNumHeavyAtoms())

        return SMILESValidationResult(
            input_smiles=smiles_clean,
            is_valid=True,
            canonical_smiles=canonical_smiles,
            molecular_formula=formula,
            molecular_weight=round(mw, 4),
            logp=round(logp, 2),
            hbd=hbd,
            hba=hba,
            tpsa=round(tpsa, 2),
            heavy_atom_count=heavy_atoms,
        )
    except Exception as e:
        return SMILESValidationResult(
            input_smiles=smiles_clean,
            is_valid=False,
            error_message=f"Error calculating properties: {e}",
        )


def generate_structure(
    smiles: str,
    output_path: Union[str, Path],
    format: str = "sdf",
    dimension: str = "3d",
    force: bool = False,
    title: Optional[str] = None,
) -> str:
    """
    Generate 2D or 3D molecular structure from a SMILES string and save to file.

    Supports:
        - Formats: 'sdf', 'mol', 'pdb'
        - Dimensions: '2d', '3d'

    Args:
        smiles: Input SMILES string.
        output_path: Destination file path.
        format: Structure format ('sdf', 'mol', 'pdb').
        dimension: Molecular coordinates ('2d' or '3d').
        force: Overwrite file if it already exists.
        title: Optional compound title or identifier to embed in the structure.

    Returns:
        Status string: 'Generated', 'Skipped (File already exists)', or error description.
    """
    if not RDKIT_AVAILABLE:
        return "RDKit is not installed; cannot generate structures"

    out_path = Path(output_path)
    if out_path.exists() and not force:
        return "Skipped (File already exists)"

    fmt = format.lower().strip()
    dim = dimension.lower().strip()

    if fmt not in ["sdf", "mol", "pdb"]:
        return f"Unsupported generation format: '{format}'. Supported formats: sdf, mol, pdb"

    if dim not in ["2d", "3d"]:
        return f"Unsupported dimension: '{dimension}'. Supported: 2d, 3d"

    if not smiles or not isinstance(smiles, str) or not smiles.strip():
        return "Empty or invalid SMILES string"

    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return f"Failed to parse SMILES '{smiles}'"

    if title:
        mol.SetProp("_Name", str(title))

    try:
        if dim == "2d":
            AllChem.Compute2DCoords(mol)
        else:
            mol = Chem.AddHs(mol)
            try:
                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                embed_res = AllChem.EmbedMolecule(mol, params)
            except AttributeError:
                embed_res = AllChem.EmbedMolecule(mol, randomSeed=42)

            if embed_res != 0:
                embed_res = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)

            if embed_res == 0:
                try:
                    AllChem.MMFFOptimizeMolecule(mol)
                except Exception:
                    try:
                        AllChem.UFFOptimizeMolecule(mol)
                    except Exception:
                        pass
            else:
                AllChem.Compute2DCoords(mol)

        out_path.parent.mkdir(parents=True, exist_ok=True)

        if fmt == "sdf":
            writer = Chem.SDWriter(str(out_path))
            writer.write(mol)
            writer.close()
        elif fmt == "mol":
            mol_block = Chem.MolToMolBlock(mol)
            out_path.write_text(mol_block, encoding="utf-8")
        elif fmt == "pdb":
            pdb_block = Chem.MolToPDBBlock(mol)
            out_path.write_text(pdb_block, encoding="utf-8")

        return "Generated"

    except Exception as e:
        return f"Generation Error: {e}"
