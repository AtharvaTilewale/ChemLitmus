"""Tests for molecular identity, library diff, and SMILES diagnosis."""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from chemlitmus import (
    DIAGNOSTIC_CATEGORIES,
    IDENTITY_LEVELS,
    compute_identity,
    describe_difference,
    diagnose_smiles,
    diff_libraries,
    group_by_identity,
    strictest_shared_level,
)
from chemlitmus.cli import app

runner = CliRunner()

L_ALA = "C[C@H](N)C(=O)O"
D_ALA = "C[C@@H](N)C(=O)O"
RAC_ALA_HCL_NA = "CC(N)C(=O)O.[Na+].[Cl-]"
RAC_ALA_NA = "CC(N)C(=O)[O-].[Na+]"
PYRIDONE_A, PYRIDONE_B = "Oc1ccccn1", "O=c1cccc[nH]1"


# --------------------------------------------------------------------------- identity

def test_identity_levels_constant():
    assert IDENTITY_LEVELS == ["exact", "parent", "tautomer", "nostereo", "skeleton", "formula"]


def test_compute_identity_keys():
    k = compute_identity(RAC_ALA_HCL_NA)
    assert k.is_valid and k.n_fragments == 3 and k.had_charge and not k.has_stereo
    assert k.exact != k.parent                      # salt stripped
    assert k.parent == "CC(N)C(=O)O"
    assert k.formula == "C3H7NO2"                   # formula of the parent, not the salt
    assert k.nostereo == k.parent                   # no stereo to remove
    for lv in IDENTITY_LEVELS:
        assert k.key(lv)
    with pytest.raises(ValueError):
        k.key("flavour")


def test_compute_identity_invalid():
    k = compute_identity("C1CC")
    assert not k.is_valid and k.error and k.exact is None


def test_identity_nesting():
    a, b = compute_identity(L_ALA), compute_identity(D_ALA)
    assert a.exact != b.exact and a.parent != b.parent
    assert a.nostereo == b.nostereo and a.skeleton == b.skeleton
    assert a.tautomer != b.tautomer                 # tautomer hash keeps stereo
    t1, t2 = compute_identity(PYRIDONE_A), compute_identity(PYRIDONE_B)
    assert t1.parent != t2.parent and t1.tautomer == t2.tautomer and t1.nostereo != t2.nostereo


def test_strictest_shared_level_and_description():
    assert strictest_shared_level(compute_identity(L_ALA), compute_identity(L_ALA)) == "exact"
    assert strictest_shared_level(compute_identity(L_ALA), compute_identity(D_ALA)) == "nostereo"
    assert describe_difference(compute_identity(L_ALA), compute_identity(D_ALA)) == "stereochemistry"
    assert describe_difference(compute_identity(RAC_ALA_HCL_NA), compute_identity(RAC_ALA_NA)) == "salt, counter-ion or charge form"
    assert describe_difference(compute_identity(PYRIDONE_A), compute_identity(PYRIDONE_B)) == "tautomer"
    assert describe_difference(compute_identity("CCO"), compute_identity("COC")).startswith("constitution")
    assert describe_difference(compute_identity("CCO"), compute_identity("c1ccccc1")) == "different compounds"
    assert strictest_shared_level(compute_identity("CCO"), compute_identity("C1CC")) is None


def test_group_by_identity():
    smiles = [L_ALA, D_ALA, RAC_ALA_HCL_NA, RAC_ALA_NA, PYRIDONE_A, PYRIDONE_B, "c1ccccc1", "CCO", "bad("]
    rep = group_by_identity(smiles, level="skeleton")
    assert rep.n_records == 9 and rep.n_valid == 8
    assert rep.n_groups == 4 and rep.n_collapsed == 4
    assert rep.n_groups_by_level["exact"] == 8 and rep.n_groups_by_level["formula"] <= 4
    sizes = sorted(g.size for g in rep.groups)
    assert sizes == [2, 4]
    ala = next(g for g in rep.groups if g.size == 4)
    assert set(ala.differs_by) == {"salt, counter-ion or charge form", "stereochemistry"}
    assert ala.distinct_exact == 4 and sorted(ala.indices) == [0, 1, 2, 3]
    pyr = next(g for g in rep.groups if g.size == 2)
    assert pyr.differs_by == ["tautomer"]
    # exact level: nothing collapses
    assert group_by_identity(smiles, level="exact").n_collapsed == 0
    with pytest.raises(ValueError):
        group_by_identity(smiles, level="flavour")


# --------------------------------------------------------------------------- diff

A = [L_ALA, RAC_ALA_HCL_NA, PYRIDONE_A, "c1ccccc1", "CCO", "CCO"]
B = [D_ALA, RAC_ALA_NA, PYRIDONE_B, "c1ccccc1", "CCO", "CCN"]


def test_diff_exact_level():
    d = diff_libraries(A, B, level="exact")
    assert (d.n_keys_a, d.n_keys_b) == (5, 6)
    assert (d.n_added, d.n_removed, d.n_unchanged, d.n_changed) == (4, 3, 2, 0)
    assert d.multiplicity_changes == 1            # CCO appears twice in A, once in B
    assert 0 < d.jaccard < 1


def test_diff_skeleton_level_classifies_changes():
    d = diff_libraries(A, B, level="skeleton")
    assert (d.n_added, d.n_removed, d.n_unchanged, d.n_changed) == (1, 0, 2, 2)
    assert d.changes_by_kind == {"salt, counter-ion or charge form; stereochemistry": 1, "tautomer": 1}
    added = [e for e in d.entries if e.status == "added"]
    assert added[0].smiles_b == ["CCN"] and added[0].smiles_a == []
    changed = {e.change: e for e in d.entries if e.status == "changed"}
    assert changed["tautomer"].smiles_a == [PYRIDONE_A] and changed["tautomer"].smiles_b == [PYRIDONE_B]
    assert not any(e.status == "unchanged" for e in d.entries)
    assert sum(1 for e in diff_libraries(A, B, level="skeleton", include_unchanged=True).entries if e.status == "unchanged") == 2


def test_diff_ordering_and_invalid_records():
    d = diff_libraries(["CCO", "bad(", "c1ccccc1"], ["CCO", "CCN"], level="parent")
    assert d.n_a == 3 and d.n_valid_a == 2
    assert [e.status for e in d.entries] == ["removed", "added"]
    with pytest.raises(ValueError):
        diff_libraries(A, B, level="flavour")


def test_diff_identical_collections():
    d = diff_libraries(A, A, level="parent")
    assert d.n_added == d.n_removed == d.n_changed == 0 and d.jaccard == 1.0


# --------------------------------------------------------------------------- diagnose

def test_diagnostic_categories_constant():
    assert DIAGNOSTIC_CATEGORIES[0] == "characters" and DIAGNOSTIC_CATEGORIES[-1] == "aromaticity"


def test_diagnose_valid():
    d = diagnose_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert d.is_valid and d.problems == [] and d.canonical_smiles and d.primary_category is None


@pytest.mark.parametrize("smiles,category,position", [
    ("C1CC(C", "parentheses", 4),
    ("C1CCC)C", "parentheses", 5),
    ("C1CCC", "rings", 1),
    ("CC1CC", "rings", 2),
    ("[Xx]C", "brackets", 0),
    ("[C@@H", "brackets", 0),
    ("CC\u2013O", "characters", 2),
    ("CC O", "characters", 2),
    ("CN(C)(C)C", "valence", 1),
    ("C(=O)(=O)(=O)C", "valence", 0),
    ("c1cccc1", "aromaticity", 0),
    ("c1cncc1", "aromaticity", 0),
    ("c1cc2ccccc2n1c", "aromaticity", 13),
])
def test_diagnose_locates_problems(smiles, category, position):
    d = diagnose_smiles(smiles)
    assert not d.is_valid
    probs = [p for p in d.problems if p.category == category]
    assert probs, [(p.category, p.message) for p in d.problems]
    assert probs[0].position == position
    assert d.primary_category in DIAGNOSTIC_CATEGORIES
    assert len(d.caret_line()) <= len(smiles) and "^" in d.caret_line()


def test_diagnose_no_duplicate_reports_for_non_ascii():
    d = diagnose_smiles("CC\u2013O")
    assert len(d.problems) == 1


def test_diagnose_valence_suggestions_and_atom_index():
    d = diagnose_smiles("CN(C)(C)C")
    p = d.problems[0]
    assert p.atom_index == 1 and "[N+]" in p.suggestion
    assert diagnose_smiles("c1cccc1").problems[0].atom_indices == [0, 1, 2, 3, 4]


@pytest.mark.parametrize("smiles,expected", [
    ("C1CC(C", "CCCC"),
    ("C1CCC)C", "CCCCC"),
    ("CC1CC", "CCCC"),
    ("CC(C)C(", "CC(C)C"),
    ("CC\u2013O", "CCO"),
    ("C C O", "CCO"),
    ("c1cncc1", "c1cc[nH]c1"),
    ("c1cc2ccccc2n1c", "Cn1ccc2ccccc21"),
    ("Cc1ccc(cc1)c", "Cc1ccc(C)cc1"),
])
def test_diagnose_repairs(smiles, expected):
    d = diagnose_smiles(smiles)
    assert d.repaired_is_valid is True and d.repaired_smiles == expected
    assert d.repairs_applied


def test_diagnose_no_repair_for_chemistry_errors_and_opt_out():
    d = diagnose_smiles("CN(C)(C)C")
    assert d.repaired_smiles is None and d.repairs_applied == []
    d2 = diagnose_smiles("C1CC(C", try_repair=False)
    assert d2.repaired_smiles is None and not d2.is_valid


def test_diagnose_empty():
    d = diagnose_smiles("   ")
    assert not d.is_valid and d.problems[0].message == "Empty SMILES"


def test_diagnose_json_roundtrip():
    d = diagnose_smiles("C1CC(C")
    data = json.loads(d.model_dump_json())
    assert data["problems"][0]["category"] in DIAGNOSTIC_CATEGORIES


# --------------------------------------------------------------------------- CLI

def _write(tmp_path: Path, name: str, lines):
    f = tmp_path / name
    f.write_text("\n".join(lines) + "\n")
    return f


def test_cli_identity_single_and_file(tmp_path: Path):
    r = runner.invoke(app, ["identity", RAC_ALA_HCL_NA])
    assert r.exit_code == 0, r.output
    assert "C3H7NO2" in r.output and "fragments=3" in r.output

    lib = _write(tmp_path, "lib.smi", [L_ALA, D_ALA, RAC_ALA_NA, PYRIDONE_A, PYRIDONE_B, "CCO"])
    out = tmp_path / "identity.csv"
    r2 = runner.invoke(app, ["identity", "--file", str(lib), "--level", "skeleton", "--output", str(out)])
    assert r2.exit_code == 0, r2.output
    assert "Identity Report" in r2.output and "stereochemistry" in r2.output
    rows = out.read_text().splitlines()
    assert len(rows) == 7 and "group_id_skeleton" in rows[0]

    assert runner.invoke(app, ["identity", "CCO", "--level", "flavour"]).exit_code == 1
    assert runner.invoke(app, ["identity"]).exit_code == 1


def test_cli_diff(tmp_path: Path):
    fa, fb = _write(tmp_path, "a.smi", A), _write(tmp_path, "b.smi", B)
    out, js = tmp_path / "diff.csv", tmp_path / "diff.json"
    r = runner.invoke(app, ["diff", str(fa), str(fb), "--level", "skeleton", "--output", str(out), "--json", str(js)])
    assert r.exit_code == 0, r.output
    assert "Library Diff" in r.output and "tautomer" in r.output
    data = json.loads(js.read_text())
    assert (data["n_added"], data["n_removed"], data["n_changed"]) == (1, 0, 2)
    assert out.read_text().count("\n") == 4  # header + 1 added + 2 changed
    assert runner.invoke(app, ["diff", str(fa), str(fb), "--level", "flavour"]).exit_code == 1


def test_cli_diagnose(tmp_path: Path):
    r = runner.invoke(app, ["diagnose", "C1CC(C"])
    assert r.exit_code == 2, r.output
    assert "parentheses" in r.output and "Candidate repair" in r.output and "CCCC" in r.output
    assert runner.invoke(app, ["diagnose", "CCO"]).exit_code == 0

    f = _write(tmp_path, "mixed.smi", ["CCO", "C1CC(C", "CN(C)(C)C", "c1cncc1"])
    out = tmp_path / "diag.csv"
    r2 = runner.invoke(app, ["diagnose", "--file", str(f), "--output", str(out)])
    assert r2.exit_code == 0, r2.output
    assert "invalid total" in r2.output
    rows = out.read_text().splitlines()
    assert len(rows) == 4  # header + 3 invalid (only-invalid default)
    r3 = runner.invoke(app, ["diagnose", "--file", str(f), "--all", "--output", str(out)])
    assert r3.exit_code == 0 and out.read_text().count("\n") == 5
    assert runner.invoke(app, ["diagnose"]).exit_code == 1
