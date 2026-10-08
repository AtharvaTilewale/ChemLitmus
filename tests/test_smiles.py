"""Tests for SMILES validation and RDKit descriptors and structure generation."""

from pathlib import Path
from chemlitmus.core.smiles import validate_smiles, generate_structure


def test_valid_smiles() -> None:
    """Test validation of valid SMILES."""
    result = validate_smiles("c1ccccc1")
    assert result.is_valid is True
    assert result.canonical_smiles == "c1ccccc1"
    assert result.molecular_formula == "C6H6"
    assert result.molecular_weight is not None
    assert result.molecular_weight > 78.0


def test_invalid_smiles() -> None:
    """Test validation of invalid SMILES."""
    result = validate_smiles("invalid_smiles_string_123")
    assert result.is_valid is False
    assert result.error_message is not None


def test_generate_structure_2d_and_3d_formats(tmp_path: Path) -> None:
    """Test generating 2D and 3D structures in SDF, MOL, and PDB formats."""
    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin

    # 2D formats
    sdf_2d = tmp_path / "aspirin_2d.sdf"
    mol_2d = tmp_path / "aspirin_2d.mol"
    pdb_2d = tmp_path / "aspirin_2d.pdb"

    assert generate_structure(smiles, sdf_2d, format="sdf", dimension="2d") == "Generated"
    assert sdf_2d.exists() and sdf_2d.stat().st_size > 0

    assert generate_structure(smiles, mol_2d, format="mol", dimension="2d") == "Generated"
    assert mol_2d.exists() and mol_2d.stat().st_size > 0

    assert generate_structure(smiles, pdb_2d, format="pdb", dimension="2d") == "Generated"
    assert pdb_2d.exists() and pdb_2d.stat().st_size > 0

    # 3D formats
    sdf_3d = tmp_path / "aspirin_3d.sdf"
    mol_3d = tmp_path / "aspirin_3d.mol"
    pdb_3d = tmp_path / "aspirin_3d.pdb"

    assert generate_structure(smiles, sdf_3d, format="sdf", dimension="3d", title="Aspirin") == "Generated"
    assert sdf_3d.exists() and sdf_3d.stat().st_size > 0

    assert generate_structure(smiles, mol_3d, format="mol", dimension="3d") == "Generated"
    assert mol_3d.exists() and mol_3d.stat().st_size > 0

    assert generate_structure(smiles, pdb_3d, format="pdb", dimension="3d") == "Generated"
    assert pdb_3d.exists() and pdb_3d.stat().st_size > 0


def test_generate_structure_resume_and_force(tmp_path: Path) -> None:
    """Test smart resume (skipping existing file) and force overwrite."""
    smiles = "CCO"
    out_file = tmp_path / "ethanol_3d.sdf"

    # First generation
    assert generate_structure(smiles, out_file, format="sdf", dimension="3d") == "Generated"
    assert out_file.exists()

    # Second generation without force should skip
    assert generate_structure(smiles, out_file, format="sdf", dimension="3d", force=False) == "Skipped (File already exists)"

    # Third generation with force should overwrite
    assert generate_structure(smiles, out_file, format="sdf", dimension="3d", force=True) == "Generated"


def test_generate_structure_invalid_inputs(tmp_path: Path) -> None:
    """Test error handling for invalid format, dimension, or SMILES."""
    out_file = tmp_path / "test.sdf"

    # Invalid SMILES
    res = generate_structure("invalid_smiles_123", out_file)
    assert "Failed to parse SMILES" in res

    # Unsupported format
    res = generate_structure("CCO", tmp_path / "test.xyz", format="xyz")
    assert "Unsupported generation format" in res

    # Unsupported dimension
    res = generate_structure("CCO", tmp_path / "test.sdf", dimension="4d")
    assert "Unsupported dimension" in res
