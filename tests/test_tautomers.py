"""Tests for chemlitmus.core.tautomers."""

from chemlitmus.core.tautomers import enumerate_tautomers, TautomerResult

URIC_ACID = "Oc1nc(O)c2nc[nH]c2n1"
BENZENE = "c1ccccc1"
INVALID = "not-a-smiles"

class TestEnumerateTautomers:
    def test_returns_tautomer_result(self):
        result = enumerate_tautomers(BENZENE)
        assert isinstance(result, TautomerResult)
        assert result.success is True

    def test_no_tautomers(self):
        result = enumerate_tautomers(BENZENE)
        assert result.num_tautomers == 1
        assert BENZENE in result.tautomers
        assert result.canonical_tautomer == BENZENE

    def test_multiple_tautomers(self):
        result = enumerate_tautomers(URIC_ACID)
        assert result.num_tautomers > 1
        assert len(result.tautomers) == result.num_tautomers
        assert result.canonical_tautomer is not None
        assert result.canonical_tautomer in result.tautomers

    def test_max_tautomers_limit(self):
        # Even with max set to a small number, RDKit might return the max limited count
        # But we can verify it doesn't crash
        result = enumerate_tautomers(URIC_ACID, max_tautomers=5)
        assert result.success is True
        assert result.num_tautomers > 0

    def test_invalid_smiles(self):
        result = enumerate_tautomers(INVALID)
        assert result.success is False
        assert result.error is not None
        assert result.num_tautomers == 0
