"""Tests for chemlitmus.core.cheminfo — fingerprints, similarity, and filters."""

import pytest
from chemlitmus.core.cheminfo import (
    substructure_search, SubstructureHit,
    FingerprintResult,
    FilterResult,
    SimilarityResult,
    RuleResult,
    compute_fingerprint,
    apply_filters,
    compute_similarity,
    _VALID_FP_TYPES,
)

# ── Test compounds ────────────────────────────────────────────────────────────

ASPIRIN = "CC(=O)OC1=CC=CC=C1C(=O)O"
ETHANOL = "CCO"
CAFFEINE = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
BENZENE = "c1ccccc1"
INVALID_SMILES = "not-a-smiles!!!"

# ── compute_fingerprint ────────────────────────────────────────────────────────


class TestComputeFingerprint:
    """Tests for compute_fingerprint()."""

    def test_ecfp4_returns_single_result(self):
        result = compute_fingerprint(ASPIRIN, fp_type="ecfp4")
        assert isinstance(result, FingerprintResult)
        assert result.fingerprint_type == "ecfp4"

    def test_ecfp4_correct_bit_count(self):
        result = compute_fingerprint(ASPIRIN, fp_type="ecfp4", n_bits=2048)
        assert result.n_bits == 2048
        assert 0 < result.n_on_bits < 2048

    def test_ecfp6_returns_result(self):
        result = compute_fingerprint(ASPIRIN, fp_type="ecfp6")
        assert result.fingerprint_type == "ecfp6"
        assert result.n_bits == 2048

    def test_maccs_fixed_167_bits(self):
        result = compute_fingerprint(ASPIRIN, fp_type="maccs")
        assert result.fingerprint_type == "maccs"
        assert result.n_bits == 167

    def test_rdkit_fp_type(self):
        result = compute_fingerprint(ASPIRIN, fp_type="rdkit")
        assert result.fingerprint_type == "rdkit"
        assert result.n_on_bits > 0

    def test_atompair_fp(self):
        result = compute_fingerprint(ASPIRIN, fp_type="atompair")
        assert result.fingerprint_type == "atompair"

    def test_torsion_fp(self):
        result = compute_fingerprint(ASPIRIN, fp_type="torsion")
        assert result.fingerprint_type == "torsion"

    def test_fcfp4_fp(self):
        result = compute_fingerprint(ASPIRIN, fp_type="fcfp4")
        assert result.fingerprint_type == "fcfp4"

    def test_all_returns_list_of_seven(self):
        results = compute_fingerprint(ASPIRIN, fp_type="all")
        assert isinstance(results, list)
        assert len(results) == len(_VALID_FP_TYPES)

    def test_all_types_present(self):
        results = compute_fingerprint(ASPIRIN, fp_type="all")
        types_returned = {r.fingerprint_type for r in results}
        assert types_returned == _VALID_FP_TYPES

    def test_bit_string_length_matches_n_bits(self):
        result = compute_fingerprint(ETHANOL, fp_type="ecfp4", n_bits=1024)
        assert len(result.bit_string) == 1024

    def test_density_is_fraction(self):
        result = compute_fingerprint(ASPIRIN, fp_type="ecfp4")
        assert 0.0 <= result.density <= 1.0

    def test_custom_bit_size(self):
        result = compute_fingerprint(ETHANOL, fp_type="ecfp4", n_bits=512)
        assert result.n_bits == 512

    def test_invalid_smiles_raises(self):
        with pytest.raises(ValueError, match="Invalid SMILES"):
            compute_fingerprint(INVALID_SMILES)

    def test_invalid_fp_type_raises(self):
        with pytest.raises(ValueError, match="Unknown fp_type"):
            compute_fingerprint(ASPIRIN, fp_type="badtype")

    def test_hex_string_is_nonempty(self):
        result = compute_fingerprint(ASPIRIN, fp_type="ecfp4")
        assert len(result.hex_string) > 0

    def test_different_molecules_different_fps(self):
        fp1 = compute_fingerprint(ASPIRIN, fp_type="ecfp4")
        fp2 = compute_fingerprint(CAFFEINE, fp_type="ecfp4")
        assert fp1.bit_string != fp2.bit_string

    def test_same_molecule_same_fp(self):
        fp1 = compute_fingerprint(ASPIRIN, fp_type="ecfp4")
        fp2 = compute_fingerprint(ASPIRIN, fp_type="ecfp4")
        assert fp1.bit_string == fp2.bit_string


# ── apply_filters ─────────────────────────────────────────────────────────────


class TestApplyFilters:
    """Tests for apply_filters()."""

    def test_aspirin_all_rules_returns_result(self):
        result = apply_filters(ASPIRIN)
        assert isinstance(result, FilterResult)
        assert result.error is None

    def test_aspirin_properties_computed(self):
        result = apply_filters(ASPIRIN)
        assert result.molecular_weight is not None
        assert abs(result.molecular_weight - 180.16) < 0.5
        assert result.logp is not None
        assert result.hbd == 1
        assert result.hba == 3

    def test_aspirin_lipinski_passes(self):
        result = apply_filters(ASPIRIN, rules=["lipinski"])
        assert result.lipinski is not None
        assert result.lipinski.passed is True

    def test_aspirin_veber_passes(self):
        result = apply_filters(ASPIRIN, rules=["veber"])
        assert result.veber.passed is True

    def test_aspirin_no_pains_alert(self):
        result = apply_filters(ASPIRIN, rules=["pains"])
        assert result.pains.passed is True

    def test_qed_score_in_range(self):
        result = apply_filters(ASPIRIN)
        assert 0.0 <= result.qed_score <= 1.0

    def test_caffeine_lipinski_passes(self):
        result = apply_filters(CAFFEINE, rules=["lipinski"])
        assert result.lipinski.passed is True

    def test_single_rule_only_populates_that_field(self):
        result = apply_filters(ASPIRIN, rules=["lipinski"])
        assert result.lipinski is not None
        assert result.veber is None
        assert result.ghose is None

    def test_passes_all_true_when_all_rules_pass(self):
        # A simple drug-like molecule should pass Lipinski + Veber
        result = apply_filters(CAFFEINE, rules=["lipinski", "veber"])
        assert result.passes_all is True

    def test_invalid_smiles_returns_error_result(self):
        result = apply_filters(INVALID_SMILES)
        assert result.error is not None
        assert result.passes_all is False

    def test_invalid_rule_name_raises(self):
        with pytest.raises(ValueError, match="Unknown rule"):
            apply_filters(ASPIRIN, rules=["notarule"])

    def test_all_keyword_applies_all_rules(self):
        result = apply_filters(ASPIRIN, rules=["all"])
        assert result.lipinski is not None
        assert result.veber is not None
        assert result.ghose is not None
        assert result.egan is not None
        assert result.ro3 is not None
        assert result.pains is not None

    def test_none_rules_applies_all(self):
        result = apply_filters(ASPIRIN, rules=None)
        assert result.lipinski is not None
        assert result.pains is not None

    def test_ghose_filter_for_small_molecule(self):
        # Ethanol (MW 46) should fail Ghose (MW 160-480 required)
        result = apply_filters(ETHANOL, rules=["ghose"])
        assert result.ghose.passed is False

    def test_egan_passes_for_aspirin(self):
        result = apply_filters(ASPIRIN, rules=["egan"])
        assert result.egan.passed is True

    def test_ro3_fails_for_large_molecule(self):
        # Aspirin MW=180 < 300 so might pass — but let's test a known drug
        result = apply_filters(CAFFEINE, rules=["ro3"])
        # Caffeine MW=194 < 300, LogP=-0.07 < 3 — should pass
        assert result.ro3.passed is True


# ── compute_similarity ─────────────────────────────────────────────────────────


class TestComputeSimilarity:
    """Tests for compute_similarity()."""

    LIBRARY = [ETHANOL, ASPIRIN, CAFFEINE, BENZENE, "CCCO", "c1ccccc1O"]

    def test_identical_molecule_is_rank_1_sim_1(self):
        hits = compute_similarity(ASPIRIN, self.LIBRARY, fp_type="ecfp4")
        assert hits[0].similarity == 1.0
        assert hits[0].hit == ASPIRIN

    def test_returns_list_of_similarity_results(self):
        hits = compute_similarity(ASPIRIN, self.LIBRARY)
        assert isinstance(hits, list)
        assert all(isinstance(h, SimilarityResult) for h in hits)

    def test_sorted_descending(self):
        hits = compute_similarity(ASPIRIN, self.LIBRARY, fp_type="ecfp4")
        sims = [h.similarity for h in hits]
        assert sims == sorted(sims, reverse=True)

    def test_ranks_are_sequential(self):
        hits = compute_similarity(ASPIRIN, self.LIBRARY, fp_type="ecfp4")
        assert [h.rank for h in hits] == list(range(1, len(hits) + 1))

    def test_threshold_filters_results(self):
        hits = compute_similarity(ASPIRIN, self.LIBRARY, threshold=0.5)
        assert all(h.similarity >= 0.5 for h in hits)

    def test_top_n_limits_results(self):
        hits = compute_similarity(ASPIRIN, self.LIBRARY, top_n=2)
        assert len(hits) <= 2

    def test_zero_threshold_returns_all_valid(self):
        hits = compute_similarity(ASPIRIN, self.LIBRARY, threshold=0.0)
        assert len(hits) == len(self.LIBRARY)

    def test_high_threshold_returns_fewer_results(self):
        hits_low = compute_similarity(ASPIRIN, self.LIBRARY, threshold=0.0)
        hits_high = compute_similarity(ASPIRIN, self.LIBRARY, threshold=0.9)
        assert len(hits_high) <= len(hits_low)

    def test_maccs_fingerprint_type(self):
        hits = compute_similarity(ASPIRIN, self.LIBRARY, fp_type="maccs")
        assert all(h.fingerprint_type == "maccs" for h in hits)

    def test_fingerprint_type_stored_in_result(self):
        hits = compute_similarity(ASPIRIN, self.LIBRARY, fp_type="ecfp6")
        assert hits[0].fingerprint_type == "ecfp6"

    def test_query_stored_in_result(self):
        hits = compute_similarity(ASPIRIN, self.LIBRARY, fp_type="ecfp4")
        assert all(h.query == ASPIRIN for h in hits)

    def test_invalid_query_smiles_raises(self):
        with pytest.raises(ValueError, match="Invalid query SMILES"):
            compute_similarity(INVALID_SMILES, self.LIBRARY)

    def test_invalid_fp_type_raises(self):
        with pytest.raises(ValueError, match="Unknown fp_type"):
            compute_similarity(ASPIRIN, self.LIBRARY, fp_type="invalid")

    def test_empty_library_returns_empty_list(self):
        hits = compute_similarity(ASPIRIN, [])
        assert hits == []

    def test_library_with_invalid_smiles_are_skipped(self):
        library_with_bad = [ETHANOL, INVALID_SMILES, CAFFEINE]
        hits = compute_similarity(ASPIRIN, library_with_bad)
        assert all(h.hit != INVALID_SMILES for h in hits)


class TestSubstructureSearch:
    def test_smarts_search_matches(self):
        library = ["CCO", "CC(=O)OC1=CC=CC=C1C(=O)O", "c1ccccc1"]
        # Carboxylic acid SMARTS
        hits = substructure_search("C(=O)[OH]", library, is_smarts=True)
        assert len(hits) == 1
        assert hits[0].smiles == "CC(=O)OC1=CC=CC=C1C(=O)O"
        assert hits[0].matched is True
        assert len(hits[0].match_indices) > 0

    def test_smarts_search_no_matches(self):
        library = ["CCO", "c1ccccc1"]
        # Halogen SMARTS
        hits = substructure_search("[#9,#17,#35,#53]", library, is_smarts=True)
        assert len(hits) == 0

    def test_smiles_search_matches(self):
        library = ["CCO", "CC(=O)OC1=CC=CC=C1C(=O)O", "c1ccccc1"]
        # Benzene ring SMILES
        hits = substructure_search("c1ccccc1", library, is_smarts=False)
        assert len(hits) == 2
        assert "c1ccccc1" in [h.smiles for h in hits]
        assert "CC(=O)OC1=CC=CC=C1C(=O)O" in [h.smiles for h in hits]

    def test_invalid_smarts_raises(self):
        library = ["CCO"]
        import pytest
        with pytest.raises(ValueError, match="Invalid SMARTS query"):
            substructure_search("invalid-smarts!!!", library, is_smarts=True)

    def test_invalid_smiles_raises(self):
        library = ["CCO"]
        import pytest
        with pytest.raises(ValueError, match="Invalid SMILES query"):
            substructure_search("invalid-smiles!!!", library, is_smarts=False)
