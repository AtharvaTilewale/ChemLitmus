"""Tests for chemlitmus.core.standardize and chemlitmus.core.iupacname."""

import pytest
from chemlitmus.core.standardize import (
    standardize_smiles,
    StandardizeResult,
    StepResult,
    VALID_STEPS,
)
from chemlitmus.core.iupacname import get_iupac_name, IUPACResult


# ── Constants ─────────────────────────────────────────────────────────────────

ASPIRIN        = "CC(=O)OC1=CC=CC=C1C(=O)O"
ETHANOL        = "CCO"
SALT_SMILES    = "[Na+].[OH-].CC(=O)[O-]"   # sodium + hydroxide + acetate
CHARGED_BA     = "O=C([O-])c1ccccc1"         # benzoate anion
PHENOL         = "Oc1ccccc1"
INVALID        = "not-a-smiles!!!"


# ═══════════════════════════════════════════════════════════════════════════════
#  standardize_smiles()
# ═══════════════════════════════════════════════════════════════════════════════

class TestStandardizeSmiles:

    # ── Return type ───────────────────────────────────────────────────────────

    def test_returns_standardize_result(self):
        result = standardize_smiles(ETHANOL)
        assert isinstance(result, StandardizeResult)

    def test_step_results_are_step_result_objects(self):
        result = standardize_smiles(ASPIRIN)
        for sr in result.step_results:
            assert isinstance(sr, StepResult)

    # ── Basic correctness ─────────────────────────────────────────────────────

    def test_already_canonical_unchanged(self):
        result = standardize_smiles(ETHANOL)
        assert result.output_smiles is not None
        assert not result.changed  # CCO is already fully canonical

    def test_salt_stripping_removes_sodium(self):
        result = standardize_smiles(SALT_SMILES)
        out = result.output_smiles or ""
        assert "[Na+]" not in out
        assert "[OH-]" not in out

    def test_salt_full_pipeline_output(self):
        result = standardize_smiles(SALT_SMILES)
        assert result.output_smiles == "CC(=O)O"
        assert result.changed is True

    def test_neutralize_removes_charge(self):
        result = standardize_smiles(CHARGED_BA, steps=["neutralize"])
        out = result.output_smiles or ""
        assert "[O-]" not in out
        assert result.changed is True

    def test_canonical_tautomer_stable(self):
        result = standardize_smiles(PHENOL, steps=["tautomer"])
        assert result.output_smiles is not None

    def test_success_property(self):
        result = standardize_smiles(ETHANOL)
        assert result.success is True

    def test_error_on_invalid_smiles(self):
        result = standardize_smiles(INVALID)
        assert result.error is not None
        assert result.output_smiles is None
        assert result.success is False

    # ── Steps configuration ───────────────────────────────────────────────────

    def test_default_applies_all_steps(self):
        result = standardize_smiles(SALT_SMILES, steps=None)
        assert result.steps_applied == VALID_STEPS

    def test_all_keyword_applies_all_steps(self):
        result = standardize_smiles(SALT_SMILES, steps=["all"])
        assert result.steps_applied == VALID_STEPS

    def test_single_step_fragment(self):
        result = standardize_smiles(SALT_SMILES, steps=["fragment"])
        assert result.steps_applied == ["fragment"]
        assert len(result.step_results) == 1

    def test_single_step_neutralize(self):
        result = standardize_smiles(CHARGED_BA, steps=["neutralize"])
        assert result.steps_applied == ["neutralize"]

    def test_ordered_steps_subset(self):
        # Even if passed out of order, they execute in canonical order
        result = standardize_smiles(SALT_SMILES, steps=["neutralize", "fragment"])
        # fragment runs before neutralize (per VALID_STEPS ordering)
        assert result.steps_applied.index("fragment") < result.steps_applied.index("neutralize")

    def test_invalid_step_raises(self):
        with pytest.raises(ValueError, match="Unknown step"):
            standardize_smiles(ETHANOL, steps=["bad_step"])

    # ── Step results ──────────────────────────────────────────────────────────

    def test_step_result_has_expected_fields(self):
        result = standardize_smiles(SALT_SMILES)
        sr = result.step_results[0]
        assert sr.step in VALID_STEPS
        assert isinstance(sr.input_smiles, str)
        assert isinstance(sr.output_smiles, str)
        assert isinstance(sr.changed, bool)

    def test_step_result_changed_flag_correct(self):
        result = standardize_smiles(SALT_SMILES)
        frag_step = next(s for s in result.step_results if s.step == "fragment")
        assert frag_step.changed is True  # salt stripped

    def test_no_steps_changed_for_clean_molecule(self):
        result = standardize_smiles(ETHANOL)
        for sr in result.step_results:
            assert sr.changed is False

    def test_input_smiles_stored(self):
        result = standardize_smiles(SALT_SMILES)
        assert result.input_smiles == SALT_SMILES


# ═══════════════════════════════════════════════════════════════════════════════
#  get_iupac_name()
# ═══════════════════════════════════════════════════════════════════════════════

class TestGetIUPACName:

    # ── Return type ───────────────────────────────────────────────────────────

    def test_returns_iupac_result(self):
        result = get_iupac_name(ASPIRIN)
        assert isinstance(result, IUPACResult)

    # ── Offline output (always populated) ────────────────────────────────────

    def test_canonical_smiles_populated(self):
        result = get_iupac_name(ASPIRIN)
        assert result.canonical_smiles is not None
        assert len(result.canonical_smiles) > 0

    def test_formula_populated(self):
        result = get_iupac_name(ASPIRIN)
        assert result.molecular_formula == "C9H8O4"

    def test_mw_populated_and_correct(self):
        result = get_iupac_name(ETHANOL)
        assert result.molecular_weight is not None
        assert abs(result.molecular_weight - 46.042) < 0.01

    def test_inchi_populated(self):
        result = get_iupac_name(ASPIRIN)
        assert result.inchi is not None
        assert result.inchi.startswith("InChI=")

    def test_inchikey_populated(self):
        result = get_iupac_name(ASPIRIN)
        assert result.inchikey == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    def test_inchikey_ethanol(self):
        result = get_iupac_name(ETHANOL)
        assert result.inchikey == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"

    def test_source_is_offline_inchi_when_no_cache(self):
        # Use a unique SMILES here to avoid picking up cache from other tests
        result = get_iupac_name("CC(=O)O", use_online=False)
        assert result.iupac_name_source == "offline_inchi"

    def test_iupac_name_none_without_online(self):
        result = get_iupac_name("CC(=O)O", use_online=False)
        assert result.iupac_name is None

    # ── Error handling ────────────────────────────────────────────────────────

    def test_invalid_smiles_returns_error_result(self):
        result = get_iupac_name(INVALID)
        assert result.error is not None
        assert result.canonical_smiles is None
        assert result.inchi is None

    def test_input_smiles_always_stored(self):
        result = get_iupac_name(ASPIRIN)
        assert result.input_smiles == ASPIRIN

    def test_invalid_smiles_input_stored(self):
        result = get_iupac_name(INVALID)
        assert result.input_smiles == INVALID

    # ── Online lookup (integration - marked) ─────────────────────────────────

    @pytest.mark.integration
    def test_online_fetches_iupac_name_ethanol(self):
        result = get_iupac_name(ETHANOL, use_online=True)
        assert result.iupac_name is not None
        assert "ethanol" in result.iupac_name.lower()
        assert result.iupac_name_source in ("pubchem", "cache")

    @pytest.mark.integration
    def test_online_fetches_aspirin_iupac(self):
        result = get_iupac_name(ASPIRIN, use_online=True)
        assert result.iupac_name is not None
        # PubChem returns "2-acetyloxybenzoic acid" or similar
        assert len(result.iupac_name) > 5
