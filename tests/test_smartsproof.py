"""Tests for static SMARTS containment proofs."""

import pytest
from typer.testing import CliRunner

from chemlitmus import audit_smarts, equivalent, prove_catalogue, satisfiable, subsumes
from chemlitmus.cli import app
from chemlitmus.core.smartsproof import Prover, parse_describe

runner = CliRunner()


@pytest.mark.parametrize("a,b,expected", [
    ("c1ccccc1", "[#6]", True), ("[#6]", "c1ccccc1", False),
    ("[CH3]C", "CC", True), ("CC", "[CH3]C", False),
    ("CCO", "CC", True), ("CC", "CCO", False),
    ("C=O", "[#6]~[#8]", True), ("[#6]~[#8]", "C=O", False),
    ("c1ccccc1[OH]", "c[OH]", True), ("c1ccccc1O", "c[OH]", False),      # O vs [OH]: H count not implied
    ("[N+](=O)[O-]", "[N+]", True), ("[C;R]", "[R]", True), ("[R]", "[C;R]", False),
    ("[nH]", "[#7;a;H1]", True), ("[#7;a;H1]", "[nH]", True),
    ("[CX4]", "[C;!$(C=O)]", None),                                       # recursive: not analysable
    ("C(=O)[OH]", "[CX3](=O)[OX2H1]", False),                             # degree of C not bounded by A: unprovable
    ("C-C", "C~C", True), ("C~C", "C-C", False), ("c:c", "c~c", True), ("C=C", "C=,:C", True),
    ("C-!@C", "CC", True), ("C!@C", "CC", False), ("CC", "C!@C", False),       # bare !@ allows any bond order
])
def test_subsumes_directional(a, b, expected):
    r = subsumes(a, b)
    if expected is None:
        assert r.status.startswith("not analysable") and r.proven is False
    else:
        assert r.proven is expected, r.status
        if expected:
            assert r.witness is not None and len(r.witness) == len(r.witness.values())   # injective


def test_witness_is_checkable():
    r = subsumes("c1ccccc1[OH]", "c[OH]")
    # B atom 1 ([OH]) must map onto A atom 6 (the hydroxyl O), B atom 0 (c) onto a ring carbon
    assert r.witness[1] == 6 and r.witness[0] in range(6)


def test_equivalent():
    assert equivalent("[CH3]", "[C;H3]").proven
    assert equivalent("[#6;a]", "c").proven
    r = equivalent("CC", "CCC")
    assert not r.proven and r.status == "no witness"
    assert equivalent("[2H]", "[H]").status.startswith("not analysable")


@pytest.mark.parametrize("smarts,status", [
    ("[C;N]", "unsatisfiable"), ("[c;!a]", "unsatisfiable"), ("[R0;x2]", "unsatisfiable"), ("[D1;x2]", "unsatisfiable"),
    ("[a;R0]", "unsatisfiable"), ("C", "satisfiable at atom and bond level"), ("[N+]#[C-]", "satisfiable at atom and bond level"),
    ("[2H]", "not analysable: isotope"), ("C.C", "not analysable: disconnected pattern"), ("[v4]", "not analysable: total valence v<n>"),
])
def test_satisfiable(smarts, status):
    assert satisfiable(smarts).status == status


def test_parse_describe_tree():
    t = parse_describe("AtomAnd\n  AtomOr\n    AtomAtomicNum 6 = val\n    AtomAtomicNum 7 = val\n  AtomInNRings -1 != val\n")
    assert t.op == "and" and t.children[0].op == "or" and len(t.children[0].children) == 2
    assert t.children[1].name == "AtomInNRings" and t.children[1].value == -1 and t.children[1].negate


def test_universe_contains_referenced_values_and_sentinels():
    pr = Prover(["[#50;D6;+3;H4]"])
    d = pr.universe.domains
    assert 50 in d["Z"] and 0 in d["Z"] and 6 in d["D"] and 9 in d["D"] and 3 in d["q"] and 99 in d["q"] and 4 in d["H"]
    assert pr.universe.size > 0


def test_prove_catalogue():
    pats = [("[Br,Cl,I][CX4]", "alkyl halide"), ("[Cl][CX4]", "chloroalkane"), ("[#6][Cl]", "any C-Cl"),
            ("[CH3]", "methyl"), ("[C;H3]", "methyl2"), ("[$(C=O)]", "recursive"), ("[C;N]", "impossible"), ("CCCCCCCCCCCCCC", "long")]
    res = prove_catalogue(pats)
    by = {p.name: p for p in res.patterns}
    assert by["chloroalkane"].status == "proven redundant" and set(by["chloroalkane"].proven_subsumed_by) == {"alkyl halide", "any C-Cl"}
    assert by["alkyl halide"].status == "no witness"                         # [Br,Cl,I] is not contained in [Cl]
    assert by["methyl"].status == "proven equivalent" and by["methyl"].proven_equivalent_to == ["methyl2"]
    assert by["recursive"].status == "not analysable" and "recursive" in by["recursive"].reason
    assert by["impossible"].status == "unsatisfiable"
    assert res.n_not_analysable == 1 and res.n_unsatisfiable == 1 and res.n_proven_redundant == 3 and res.n_proven_equivalent == 2
    assert res.not_analysable_reasons == {"recursive SMARTS $(...)": 1}
    small = prove_catalogue(pats, max_container_atoms=1)
    assert small.n_pairs_budget_exceeded > 0 and small.n_proven_redundant <= res.n_proven_redundant


def test_audit_proof_check_is_opt_in():
    pats = ["[Cl][CX4]", "[#6][Cl]", "[CH3]", "[C;H3]"]
    res = audit_smarts(pats, checks=["compile", "proof"])
    assert res.checks_run == ["compile", "proof"]
    assert res.patterns[0].proof_status == "proven redundant" and res.patterns[0].proven_subsumed_by == [1]
    assert res.patterns[2].proof_status == "proven equivalent" and res.patterns[2].proven_equivalent_to == [3]
    assert "proven-redundant" in res.patterns[0].flags and res.n_proven_redundant == 3
    assert res.to_rows()[0]["proof_status"] == "proven redundant"
    plain = audit_smarts(pats, checks=["compile"])
    assert plain.patterns[0].proof_status is None and "proven-redundant" not in plain.patterns[0].flags


def test_cli_smartsproof(tmp_path):
    r = runner.invoke(app, ["smartsproof", "--subsumes", "c1ccccc1[OH]", "c[OH]"])
    assert r.exit_code == 0 and "proven" in r.output and "Witness" in r.output
    r = runner.invoke(app, ["smartsproof", "--subsumes", "CC", "CCO"])
    assert r.exit_code == 1 and "no witness" in r.output and "not a refutation" in r.output
    r = runner.invoke(app, ["smartsproof", "--equivalent", "[CH3]", "[C;H3]"])
    assert r.exit_code == 0 and "proven" in r.output
    r = runner.invoke(app, ["smartsproof", "--satisfiable", "[R0;x2]"])
    assert r.exit_code == 1 and "unsatisfiable" in r.output
    assert runner.invoke(app, ["smartsproof", "--subsumes", "C((", "C"]).exit_code == 1
    assert runner.invoke(app, ["smartsproof"]).exit_code == 1
    f = tmp_path / "p.csv"; f.write_text("name,smarts\nchloro,[Cl][CX4]\nany,[#6][Cl]\nrec,[$(C=O)]\n")
    out, js = tmp_path / "p_out.csv", tmp_path / "p.json"
    r = runner.invoke(app, ["smartsproof", str(f), "-o", str(out), "--json", str(js)])
    assert r.exit_code == 0 and "Proven redundant" in r.output and "chloro" in r.output
    rows = out.read_text().splitlines(); assert len(rows) == 4 and rows[1].split(",")[3] == "proven redundant"
    r = runner.invoke(app, ["smartsaudit", str(f), "--checks", "compile,proof", "--max-molecules", "50"])
    assert r.exit_code == 0 and "Proven redundant" in r.output
    assert runner.invoke(app, ["smartsaudit", str(f), "--checks", "nope"]).exit_code == 1
