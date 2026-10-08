"""Tests for the SMARTS auditing module and the --prep option on filter/substructure."""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from rdkit import Chem

from chemlitmus import (
    AUDIT_CHECKS,
    PREPARATIONS,
    SmartsAuditResult,
    SmartsExplanation,
    apply_filters,
    audit_smarts,
    explain_smarts,
    load_patterns,
    load_reference_library,
    prepare_molecule,
    substructure_search,
)
from chemlitmus.cli import app

runner = CliRunner()

# A small, fully controlled reference library.
LIBRARY_SMILES = [
    "c1ccccc1", "Cc1ccccc1", "c1ccccc1O", "c1ccncc1", "c1ccc2ccccc2c1",          # aromatics
    "CCO", "CC(C)O", "CCCO", "OCCO", "CC(O)CO",                                      # alcohols
    "CC(=O)O", "CC(=O)OC", "CC(=O)N", "CC(=O)NC", "O=C(O)c1ccccc1",                  # carbonyls
    "CCN", "CCNC", "CCN(C)C", "NCCO", "CNC(=O)C",                                    # amines/amides
    "CCCl", "CCBr", "CCF", "ClCCCl", "FC(F)F",                                       # halides
    "C1CCCCC1", "C1CCOC1", "C1CCNC1", "CCCCCC", "CC(C)(C)C",                         # aliphatics
]
LIBRARY = [Chem.MolFromSmiles(s) for s in LIBRARY_SMILES]
N = len(LIBRARY)

PATTERNS = [
    ("c1ccccc1",   "benzene",        "demo"),   # 0: fires; subsumed by [#6]
    ("[#6]",       "any carbon",     "demo"),   # 1: over-broad
    ("c1ccccc1",   "benzene again",  "demo"),   # 2: exact duplicate of 0
    ("[OX2H]",     "hydroxyl",       "demo"),   # 3: fires
    ("[OH]",       "hydroxyl alt",   "demo"),   # 4: library-equivalent to 3
    ("[#118]",     "oganesson",      "demo"),   # 5: dead, never-matching atom
    ("ClCCCCCCCl", "dichloroheptane","demo"),   # 6: dead, rare combination
    ("c[H]",       "aromatic C-H",   "demo"),   # 7: needs explicit H; fires only with explicit-h
    ("C1=CC=CC=C1","kekule benzene", "demo"),   # 8: fires only with kekule
    ("[C",         "broken",         "demo"),   # 9: unparseable
    ("[$(C=O)]N",  "amide-ish",      "demo"),   # 10: recursive
]


@pytest.fixture(scope="module")
def result() -> SmartsAuditResult:
    return audit_smarts(PATTERNS, library=LIBRARY, library_source="inline", breadth_threshold=0.5)


def by_name(res: SmartsAuditResult, name: str):
    return next(p for p in res.patterns if p.name == name)


# --------------------------------------------------------------------------- audit

def test_constants():
    assert AUDIT_CHECKS == ["compile", "breadth", "dead", "redundancy", "sensitivity"]
    assert PREPARATIONS == ["implicit-h", "explicit-h", "kekule"]


def test_shape_and_bookkeeping(result):
    assert result.n_patterns == len(PATTERNS)
    assert result.n_molecules == N
    assert result.library_source == "inline"
    assert result.checks_run == AUDIT_CHECKS
    assert result.error is None
    assert result.elapsed_seconds >= 0


def test_compile_flags(result):
    broken = by_name(result, "broken")
    assert broken.parses is False and broken.parse_error
    assert broken.n_hits is None and "unparseable" in broken.flags
    assert by_name(result, "aromatic C-H").requires_explicit_h is True
    assert by_name(result, "benzene").requires_explicit_h is False
    assert by_name(result, "amide-ish").has_recursive_smarts is True
    assert by_name(result, "hydroxyl").has_recursive_smarts is False
    assert result.n_unparseable == 1


def test_hit_counts_match_rdkit(result):
    q = Chem.MolFromSmarts("[OX2H]")
    expected = sum(1 for m in LIBRARY if m.HasSubstructMatch(q))
    assert by_name(result, "hydroxyl").n_hits == expected
    assert by_name(result, "hydroxyl").hit_fraction == pytest.approx(expected / N)


def test_breadth(result):
    assert by_name(result, "any carbon").over_broad is True
    assert by_name(result, "hydroxyl").over_broad is False
    assert result.n_over_broad >= 1


def test_dead_triage(result):
    og = by_name(result, "oganesson")
    assert og.never_fires and og.dead_verdict == "never-matching atom" and og.never_matching_atoms == [0]
    rare = by_name(result, "dichloroheptane")
    assert rare.never_fires and rare.dead_verdict == "rare combination" and rare.never_matching_atoms == []
    assert by_name(result, "benzene").never_fires is False
    assert result.n_dead_never_matching_atom == 1


def test_dead_but_alive_under_other_preparation(result):
    ch = by_name(result, "aromatic C-H")
    assert ch.never_fires and ch.dead_verdict == "fires only with explicit-h"
    assert ch.hits_by_preparation["implicit-h"] == 0 and ch.hits_by_preparation["explicit-h"] > 0
    kek = by_name(result, "kekule benzene")
    assert kek.dead_verdict == "fires only with kekule"
    assert kek.hits_by_preparation["kekule"] > 0


def test_redundancy(result):
    assert by_name(result, "benzene again").duplicate_of == 0
    assert by_name(result, "benzene").duplicate_of is None
    h1, h2 = by_name(result, "hydroxyl"), by_name(result, "hydroxyl alt")
    assert h2.index in h1.equivalent_to and h1.index in h2.equivalent_to
    assert by_name(result, "benzene").subsumed_by == 1          # every benzene contains a carbon
    assert by_name(result, "any carbon").subsumed_by is None    # nothing contains "any carbon"
    assert result.n_duplicates == 1


def test_sensitivity_summary(result):
    s = result.sensitivity
    assert s is not None and s.preparations == PREPARATIONS and s.n_molecules == N
    assert set(s.compounds_flagged) == set(PREPARATIONS)
    assert set(s.verdict_flips) == {"explicit-h", "kekule"}
    assert by_name(result, "aromatic C-H").preparation_sensitive is True
    assert by_name(result, "hydroxyl").preparation_sensitive is False
    assert s.n_sensitive_patterns >= 2


def test_flags_and_clean(result):
    assert set(by_name(result, "hydroxyl").flags) == {"equivalent", "subsumed"}  # subsumed by [#6]
    assert result.n_clean == sum(1 for p in result.patterns if not p.flags)
    rows = result.to_rows()
    assert len(rows) == len(PATTERNS)
    assert {"smarts", "flags", "hits_implicit-h", "hits_explicit-h", "hits_kekule"} <= set(rows[0])


def test_plain_string_input_and_check_subset():
    res = audit_smarts(["c1ccccc1", "[OH]"], library=LIBRARY, checks=["breadth"])
    assert res.checks_run == ["compile", "breadth"]
    assert res.sensitivity is None
    assert res.patterns[0].name is None and res.patterns[0].n_hits > 0
    assert res.patterns[0].dead_verdict is None and res.patterns[0].subsumed_by is None


def test_invalid_check_and_preparation():
    with pytest.raises(ValueError):
        audit_smarts(["C"], library=LIBRARY, checks=["nonsense"])
    with pytest.raises(ValueError):
        audit_smarts(["C"], library=LIBRARY, preparations=["gas-phase"])
    with pytest.raises(ValueError):
        prepare_molecule(LIBRARY[0], "gas-phase")


def test_json_roundtrip(result):
    data = json.loads(result.model_dump_json())
    assert data["n_patterns"] == len(PATTERNS)
    assert SmartsAuditResult.model_validate(data).n_dead == result.n_dead


# --------------------------------------------------------------------------- explain

def test_explain_fires():
    e = explain_smarts("[OX2H]", library=LIBRARY)
    assert isinstance(e, SmartsExplanation) and e.parses
    assert e.verdict == "fires" and e.n_query_atoms == 1
    assert e.hits_by_preparation["implicit-h"] > 0
    assert e.example_matches and all(Chem.MolFromSmiles(s) for s in e.example_matches)


def test_explain_dead_atom_and_explicit_h():
    e = explain_smarts("[#118]C", library=LIBRARY)
    assert e.verdict == "dead: never-matching atom" and e.never_matching_atoms == [0]
    assert e.atoms[0].n_matching_molecules == 0 and e.atoms[1].n_matching_molecules > 0
    e2 = explain_smarts("c[H]", library=LIBRARY)
    assert e2.requires_explicit_h and e2.verdict == "fires only with explicit-h"
    assert e2.atoms[1].is_hydrogen is True


def test_explain_unparseable():
    e = explain_smarts("[C", library=LIBRARY)
    assert e.parses is False and e.verdict == "unparseable" and e.parse_error


# --------------------------------------------------------------------------- loading

def test_load_patterns_csv_and_text(tmp_path: Path):
    csv = tmp_path / "alerts.csv"
    csv.write_text("rule_id,description,smarts,rule_set_name\n1,phenol,c[OH],demo\n2,,[NX3],demo\n3,blank,,demo\n")
    rows = load_patterns(csv)
    assert rows == [("c[OH]", "phenol", "demo"), ("[NX3]", None, "demo")]

    txt = tmp_path / "alerts.smarts"
    txt.write_text("# comment\nc[OH] phenol group\n\n[NX3]\n")
    assert load_patterns(txt) == [("c[OH]", "phenol group", None), ("[NX3]", None, None)]

    bad = tmp_path / "bad.csv"
    bad.write_text("a,b\n1,2\n")
    with pytest.raises(ValueError):
        load_patterns(bad)
    with pytest.raises(FileNotFoundError):
        load_patterns(tmp_path / "missing.csv")


def test_bundled_reference_library():
    mols, source = load_reference_library(max_molecules=50)
    assert len(mols) == 50 and "bundled" in source
    mols_all, _ = load_reference_library()
    assert len(mols_all) > 9000


def test_user_reference_library(tmp_path: Path):
    smi = tmp_path / "lib.smi"
    smi.write_text("\n".join(LIBRARY_SMILES))
    mols, source = load_reference_library(smi)
    assert len(mols) == N and source == str(smi)


# --------------------------------------------------------------------------- --prep on other commands

def test_prepare_molecule_states():
    m = Chem.MolFromSmiles("c1ccccc1")
    assert prepare_molecule(m, "implicit-h").GetNumAtoms() == 6
    assert prepare_molecule(m, "explicit-h").GetNumAtoms() == 12
    k = prepare_molecule(m, "kekule")
    assert not any(a.GetIsAromatic() for a in k.GetAtoms())


def test_substructure_search_preparation():
    assert substructure_search("c[H]", LIBRARY_SMILES) == []
    hits = substructure_search("c[H]", LIBRARY_SMILES, preparation="explicit-h")
    assert hits and all(h.preparation == "explicit-h" for h in hits)
    assert substructure_search("C1=CC=CC=C1", LIBRARY_SMILES) == []
    assert substructure_search("C1=CC=CC=C1", LIBRARY_SMILES, preparation="kekule")
    with pytest.raises(ValueError):
        substructure_search("C", LIBRARY_SMILES, preparation="gas-phase")


def test_apply_filters_records_preparation():
    r = apply_filters("CC(=O)Oc1ccccc1C(=O)O", rules=["pains"], preparation="explicit-h")
    assert r.preparation == "explicit-h" and r.pains is not None
    assert apply_filters("CCO", rules=["lipinski"]).preparation == "implicit-h"
    with pytest.raises(ValueError):
        apply_filters("CCO", rules=["pains"], preparation="gas-phase")


# --------------------------------------------------------------------------- CLI

def test_cli_smartsaudit_file(tmp_path: Path):
    pat = tmp_path / "alerts.smarts"
    pat.write_text("c1ccccc1 benzene\n[OH] hydroxyl\n[#118] oganesson\nc[H] aromatic-CH\n")
    lib = tmp_path / "lib.smi"
    lib.write_text("\n".join(LIBRARY_SMILES))
    out_csv, out_json = tmp_path / "audit.csv", tmp_path / "audit.json"
    r = runner.invoke(app, ["smartsaudit", str(pat), "--library", str(lib), "--output", str(out_csv), "--json", str(out_json)])
    assert r.exit_code == 0, r.output
    assert "SMARTS Audit" in r.output and "Reproducibility" in r.output
    assert out_csv.exists() and out_json.exists()
    data = json.loads(out_json.read_text())
    assert data["n_patterns"] == 4 and data["n_molecules"] == N
    names = {p["name"]: p for p in data["patterns"]}
    assert names["oganesson"]["dead_verdict"] == "never-matching atom"
    assert names["aromatic-CH"]["dead_verdict"] == "fires only with explicit-h"
    assert out_csv.read_text().count("\n") == 5  # header + 4 rows


def test_cli_smartsaudit_explain(tmp_path: Path):
    lib = tmp_path / "lib.smi"
    lib.write_text("\n".join(LIBRARY_SMILES))
    r = runner.invoke(app, ["smartsaudit", "--explain", "[#8]-[#1]", "--library", str(lib)])
    assert r.exit_code == 0, r.output
    assert "[#8]" in r.output and "[#1]" in r.output        # markup must not swallow bracket atoms
    assert "fires only with explicit-h" in r.output


def test_cli_smartsaudit_requires_input():
    r = runner.invoke(app, ["smartsaudit"])
    assert r.exit_code == 1 and "Provide a pattern file" in r.output


def test_cli_prep_option_on_substructure_and_filter(tmp_path: Path):
    lib = tmp_path / "lib.smi"
    lib.write_text("\n".join(LIBRARY_SMILES))
    out = tmp_path / "hits.csv"
    r = runner.invoke(app, ["substructure", "c[H]", "--file", str(lib), "--prep", "explicit-h", "--output", str(out)])
    assert r.exit_code == 0, r.output
    assert "preparation" in out.read_text().splitlines()[0]
    assert all(line.endswith("explicit-h") for line in out.read_text().splitlines()[1:])

    r2 = runner.invoke(app, ["filter", "c1ccccc1O", "--rules", "pains", "--prep", "kekule"])
    assert r2.exit_code == 0, r2.output
    assert "kekule" in r2.output
