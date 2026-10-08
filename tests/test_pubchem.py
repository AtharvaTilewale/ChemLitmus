"""Integration tests for PubChem REST API client."""

import pytest
from chemlitmus.core.pubchem import PubChemClient

@pytest.mark.integration
def test_pubchem_lookup_valid():
    """Test live PubChem lookup for a known compound."""
    client = PubChemClient()
    compound = client.lookup("c1ccccc1", search_type="smiles")
    
    assert compound is not None
    assert compound.cid == 241
    assert compound.molecular_formula == "C6H6"
    assert compound.iupac_name.lower() == "benzene"
    assert compound.canonical_smiles is not None

@pytest.mark.integration
def test_pubchem_lookup_invalid():
    """Test live PubChem lookup for an invalid/non-existent compound."""
    client = PubChemClient()
    compound = client.lookup("InvalidSMILES12345", search_type="smiles")
    assert compound is None

@pytest.mark.integration
def test_pubchem_lookup_by_cid():
    """Test live PubChem lookup by explicit CID."""
    client = PubChemClient()
    compound = client.lookup(2244, search_type="cid")
    
    assert compound is not None
    assert compound.cid == 2244
    assert compound.iupac_name.lower() == "2-acetyloxybenzoic acid"
    assert compound.canonical_smiles == "CC(=O)OC1=CC=CC=C1C(=O)O"

@pytest.mark.integration
def test_pubchem_lookup_by_name():
    """Test live PubChem lookup by compound name."""
    client = PubChemClient()
    compound = client.lookup_by_name("aspirin")
    
    assert compound is not None
    assert compound.cid == 2244
    assert compound.canonical_smiles == "CC(=O)OC1=CC=CC=C1C(=O)O"
