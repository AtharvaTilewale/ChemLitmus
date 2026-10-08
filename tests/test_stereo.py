"""Tests for chemlitmus.core.stereo"""
from chemlitmus.core.stereo import analyze_stereochemistry

def test_analyze_stereocenters():
    # (R)-2-butanol
    result = analyze_stereochemistry("C[C@H](O)CC")
    assert result.success is True
    assert len(result.chiral_centers) == 1
    assert result.has_unassigned is False
    assert result.chiral_centers[0]["config"] == "S" # RDKit parses this specific SMILES as S depending on canon

def test_analyze_unassigned():
    # 2-butanol without stereo marks
    result = analyze_stereochemistry("CC(O)CC")
    assert result.success is True
    assert len(result.chiral_centers) == 1
    assert result.has_unassigned is True
    assert result.chiral_centers[0]["config"] == "?"

def test_no_stereocenters():
    result = analyze_stereochemistry("CCO")
    assert result.success is True
    assert len(result.chiral_centers) == 0
    assert result.has_unassigned is False

def test_invalid_smiles():
    result = analyze_stereochemistry("invalid")
    assert result.success is False
    assert result.error is not None
