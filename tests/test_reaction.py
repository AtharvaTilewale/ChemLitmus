"""Tests for chemlitmus.core.reaction"""
from chemlitmus.core.reaction import validate_reaction

def test_valid_reaction():
    # Esterification
    result = validate_reaction("CC(=O)O.OCC>>CC(=O)OCC.O")
    assert result.is_valid is True
    assert result.num_reactants == 2
    assert result.num_products == 2
    assert result.num_agents == 0

def test_reaction_with_agents():
    # Reaction with agent/catalyst
    result = validate_reaction("C.O>[H+]>CO")
    assert result.is_valid is True
    assert result.num_reactants == 2
    assert result.num_products == 1
    assert result.num_agents == 1

def test_invalid_reaction():
    result = validate_reaction("invalid>>reaction")
    assert result.is_valid is False
    assert result.error is not None
