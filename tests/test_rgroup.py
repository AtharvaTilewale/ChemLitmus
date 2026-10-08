"""Tests for chemlitmus.core.rgroup"""
from chemlitmus.core.rgroup import rgroup_decomposition, RGroupResult

def test_rgroup_match():
    core = "c1ccccc1"
    smiles = ["Cc1ccccc1", "c1ccccc1F"]
    results = rgroup_decomposition(core, smiles)
    assert len(results) == 2
    assert results[0].is_matched is True
    assert results[1].is_matched is True
    
    # Check that decomposition has keys Core and R1
    assert "Core" in results[0].decomposition
    assert "R1" in results[0].decomposition

def test_rgroup_no_match():
    core = "c1ccccc1"
    smiles = ["CCO"]
    results = rgroup_decomposition(core, smiles)
    assert len(results) == 1
    assert results[0].is_matched is False
    assert results[0].error == "Molecule does not match core."

def test_invalid_core():
    results = rgroup_decomposition("invalid", ["CCO"])
    assert len(results) == 1
    assert results[0].is_matched is False
    assert results[0].error == "Invalid core SMARTS string."
