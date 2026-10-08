"""Tests for chemlitmus.core.scaffold"""
from chemlitmus.core.scaffold import extract_scaffold, ScaffoldResult

def test_extract_scaffold_aspirin():
    result = extract_scaffold("CC(=O)OC1=CC=CC=C1C(=O)O")
    assert result.success is True
    assert result.scaffold_smiles == "c1ccccc1"

def test_extract_scaffold_no_rings():
    # Ethanol has no rings, scaffold is empty string
    result = extract_scaffold("CCO")
    assert result.success is True
    assert result.scaffold_smiles == ""

def test_invalid_smiles():
    result = extract_scaffold("invalid")
    assert result.success is False
    assert result.error is not None
