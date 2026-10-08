"""Tests for chemlitmus.core.atommap"""
from chemlitmus.core.atommap import map_atoms

def test_map_atoms_valid():
    result = map_atoms("CCO")
    assert result.success is True
    assert result.mapped_smiles is not None
    assert ":1" in result.mapped_smiles
    assert ":2" in result.mapped_smiles
    assert ":3" in result.mapped_smiles

def test_map_atoms_invalid():
    result = map_atoms("invalid")
    assert result.success is False
    assert result.error is not None
