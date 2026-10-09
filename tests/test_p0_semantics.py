"""Regression tests for the P0 correctness fixes: fingerprint semantics, shared SMILES
validation, and explicit failure statuses in the SMARTS audit."""

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator as G

from chemlitmus.core import smartsaudit as sa
from chemlitmus.core.cheminfo import FP_DESCRIPTIONS, compute_fingerprint, compute_similarity
from chemlitmus.core.scaffold import extract_scaffold
from chemlitmus.core.smiles import SmilesParseError, is_cxsmiles, mol_from_smiles, split_smiles_field
from chemlitmus.core.standardize import standardize_smiles
from chemlitmus.core.stereo import analyze_stereochemistry
from chemlitmus.core.tautomers import enumerate_tautomers

# ------------------------------------------------------------------ 4.1 fingerprint semantics

@pytest.mark.parametrize("smi", ["c1ccccc1C(=O)O", "CCN(CC)CC", "O=C(O)C[C@H](N)C(=O)O", "Clc1ccc(cc1)C(=O)Nc1ccncc1"])
def test_fcfp4_uses_feature_invariants(smi):
    m = Chem.MolFromSmiles(smi)
    ref_f = G.GetMorganGenerator(radius=2, fpSize=2048, atomInvariantsGenerator=G.GetMorganFeatureAtomInvGen()).GetFingerprint(m).ToBitString()
    ref_e = G.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(m).ToBitString()
    f = compute_fingerprint(smi, "fcfp4"); e = compute_fingerprint(smi, "ecfp4")
    assert f.bit_string == ref_f and e.bit_string == ref_e
    assert f.atom_invariants == "feature" and e.atom_invariants == "connectivity"
    assert f.radius == e.radius == 2 and f.algorithm == "morgan" and f.representation == "bit" and f.chirality is False


def test_fingerprint_descriptions_cover_every_type():
    for t, meta in FP_DESCRIPTIONS.items():
        r = compute_fingerprint("CCO", t)
        assert r.algorithm == meta["algorithm"] and r.radius == meta.get("radius")
    # FCFP and ECFP agree on some molecules and differ on most; assert only the aggregate
    smis = ["c1ccccc1O", "CCO", "c1ccncc1", "CC(=O)Nc1ccc(O)cc1", "OC(=O)CCC(=O)O", "C1CCNCC1"]
    diff = sum(compute_fingerprint(s, "fcfp4").bit_string != compute_fingerprint(s, "ecfp4").bit_string for s in smis)
    assert diff >= len(smis) // 2


# ------------------------------------------------------------------ 4.2 shared validation

def test_mol_from_smiles_rejects_internal_whitespace():
    with pytest.raises(SmilesParseError):
        mol_from_smiles("CC O")
    with pytest.raises(SmilesParseError):
        mol_from_smiles("   ")
    assert mol_from_smiles("C((") is None
    assert Chem.MolToSmiles(mol_from_smiles(" CCO ")) == "CCO"


def test_cxsmiles_is_the_sanctioned_whitespace():
    assert is_cxsmiles("CCO |$;;R1$|") and not is_cxsmiles("CC O") and not is_cxsmiles("CCO |a| |b|")
    m = mol_from_smiles("C[C@H](N)O |&1:1|")
    assert m is not None and m.GetNumAtoms() == 4


def test_split_smiles_field():
    assert split_smiles_field("CCO ethanol") == ("CCO", "ethanol")
    assert split_smiles_field("CCO") == ("CCO", None)
    assert split_smiles_field("CCO  my name") == ("CCO", "my name")


@pytest.mark.parametrize("fn,attr", [
    (standardize_smiles, "error"), (enumerate_tautomers, "error"), (analyze_stereochemistry, "error"), (extract_scaffold, "error"),
])
def test_structure_consumers_do_not_truncate(fn, attr):
    """'CC O' must be an error everywhere — never silently ethane."""
    res = fn("CC O")
    assert getattr(res, attr), f"{fn.__name__} accepted 'CC O'"
    ok = fn("CCO")
    assert not getattr(ok, attr)


def test_fingerprint_and_similarity_reject_whitespace():
    with pytest.raises(ValueError):
        compute_fingerprint("CC O", "ecfp4")
    with pytest.raises(ValueError):
        compute_similarity("CC O", ["CCO"])


# ------------------------------------------------------------------ 4.3 failures are not negatives

LIB = [Chem.MolFromSmiles(s) for s in ["CCO", "c1ccccc1", "CC(=O)O", "c1ccncc1", "CCN"]]


def test_preparation_failure_is_counted_not_hidden(monkeypatch):
    real = sa.prepare_molecule_status

    def flaky(mol, prep):
        if prep == "kekule" and mol.GetNumAtoms() == 6 and mol.GetAtomWithIdx(0).GetIsAromatic():
            return None, "kekulization failed: KekulizeException"
        return real(mol, prep)
    monkeypatch.setattr(sa, "prepare_molecule_status", flaky)
    lib = sa.PreparedLibrary(LIB, "kekule")
    assert lib.n_evaluated == 3 and lib.n_failed == 2 and lib.failures == {"kekulization failed: KekulizeException": 2}
    v, err = lib.match(Chem.MolFromSmarts("[#6]"))
    assert err is None and v.tolist() == [True, False, True, False, True]        # unevaluated molecules are not hits
    res = sa.audit_smarts(["c1ccccc1", "[OH]"], library=LIB, checks=["compile", "breadth", "dead", "sensitivity"])
    ps = res.preparation_status["kekule"]
    assert ps.n_failed == 2 and ps.n_evaluated == 3 and res.preparation_status["implicit-h"].n_failed == 0
    # strict preparation raises; lenient returns the original (documented legacy behaviour)
    with pytest.raises(sa.PreparationError):
        sa.prepare_molecule(LIB[1], "kekule", strict=True)
    assert sa.prepare_molecule(LIB[1], "kekule").GetNumAtoms() == 6


def test_match_error_is_reported_not_zero(monkeypatch):
    real_match = sa.PreparedLibrary.match

    def boom(self, query):
        if query.GetNumAtoms() == 1 and query.GetAtomWithIdx(0).GetAtomicNum() == 8:
            return np.zeros(len(self.evaluated), bool), "RuntimeError: simulated matcher failure"
        return real_match(self, query)
    monkeypatch.setattr(sa.PreparedLibrary, "match", boom)
    res = sa.audit_smarts(["[#8]", "[#6]"], library=LIB, checks=["compile", "breadth", "dead", "sensitivity"])
    bad, good = res.patterns
    assert bad.match_error and bad.n_hits is None and bad.never_fires is False and bad.dead_verdict is None
    assert "match-failed" in bad.flags and "dead:never-matching-atom" not in bad.flags
    assert good.n_hits == 5 and res.n_match_failures == 1
    assert res.to_rows()[0]["match_error"].startswith("implicit-h: RuntimeError")
    exp = sa.explain_smarts("[#8]", library=LIB)
    assert exp.match_errors and exp.hits_by_preparation["implicit-h"] == 0
