"""The documented commands run on the files in examples/ and find what examples/README.md says they find."""

import csv
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from chemlitmus.cli import app

EX = Path(__file__).resolve().parents[1] / "examples"
runner = CliRunner()
BIO = ["--endpoint-column", "standard_value", "--units-column", "standard_units", "--relation-column", "standard_relation"]


def _run(*args):
    r = runner.invoke(app, [str(a) for a in args])
    assert r.exit_code in (0, 3), r.output
    return r


def test_example_files_present_and_sized():
    sizes = {"bioactivity.csv": 33, "library.csv": 24, "library_v2.csv": 21, "demo.csv": 9, "demo_v2.csv": 5, "alerts.csv": 18, "alerts_v2.csv": 13}
    for name, n in sizes.items():
        assert len(list(csv.DictReader(open(EX / name)))) == n, name
    assert len((EX / "generated.smi").read_text().splitlines()) == 20
    assert len((EX / "train.smi").read_text().splitlines()) == 15
    assert len((EX / "alerts.smarts").read_text().splitlines()) == 5
    assert (EX / "README.md").exists() and (EX / "policy.json").exists()


def test_audit_bioactivity(tmp_path):
    _run("audit", EX / "bioactivity.csv", *BIO, "--split-column", "split", "--date-column", "year", "-o", tmp_path)
    a = json.loads((tmp_path / "audit.json").read_text())
    s = a["summary"]
    assert (s["n_total"], s["n_ok"], s["n_empty"], s["n_invalid"]) == (33, 30, 1, 2) and s["processing_complete"]
    codes = s["issues_by_code"]
    for c in ("SPLIT_OVERLAP", "LABEL_CONFLICT", "DUP_PARENT", "DUP_EXACT", "DUP_NOSTEREO", "FRAG_MULTIPLE", "ELEMENT_UNUSUAL", "ISOTOPE",
              "STEREO_PARTIAL", "CHARGE_NET", "ALERT_MATCH", "UNITS_MISSING", "SPLIT_MISSING", "PARSE_WHITESPACE", "PARSE_INVALID", "PARSE_EMPTY"):
        assert codes.get(c, 0) >= 1, c
    # gates: two invalid records (> 2 %) and a parent-level overlap -> exit 3
    r = runner.invoke(app, ["audit", str(EX / "bioactivity.csv"), "--split-column", "split", "--fail-on-severity", "error", "--max-invalid-fraction", "0.02", "--no-split-overlap", "-o", str(tmp_path / "g")])
    assert r.exit_code == 3
    # the example policy is loadable and changes the identity level
    r = _run("audit", EX / "bioactivity.csv", "--config", EX / "policy.json", "-o", tmp_path / "p")
    assert json.loads((tmp_path / "p" / "policy.json").read_text())["identity_level"] == "nostereo"


def test_leakage_conflicts_cliffs(tmp_path):
    _run("leakage", EX / "bioactivity.csv", "--split-column", "split", "--json", tmp_path / "l.json")
    pair = json.loads((tmp_path / "l.json").read_text())["pairs"][0]
    by = {o["level"]: o["n_eval_records"] for o in pair["overlap"]}
    assert by["exact"] == 1 and by["parent"] == 2 and by["tautomer"] == 3 and by["nostereo"] == 3 and by["skeleton"] == 4
    _run("conflicts", EX / "bioactivity.csv", *BIO, "--context", "assay_id", "--source-column", "document_id", "--json", tmp_path / "c.json")
    c = json.loads((tmp_path / "c.json").read_text())
    assert c["n_groups_compared"] == 2 and c["n_conflicts"] == 1 and c["n_missing_units"] == 1
    g = c["groups"][0]
    assert {m["source_id"] for m in g["measurements"]} == {"CPD-017", "CPD-018"} and g["spread"] == pytest.approx(3.2041, abs=1e-3)
    assert g["replicates"]["basis"] == "source field" and g["replicates"]["n_sources"] == 2
    _run("cliffs", EX / "bioactivity.csv", *BIO, "--context", "assay_id", "--json", tmp_path / "k.json", "-o", tmp_path / "pairs.csv", "--outliers", tmp_path / "out.csv")
    k = json.loads((tmp_path / "k.json").read_text())
    assert k["n_internal_conflict_excluded"] == 1 and k["n_undetermined"] >= 1 and k["n_cliffs"] >= 10
    assert [o["source_ids"] for o in k["outliers"]] == [["CPD-007"]]
    o = k["outliers"][0]
    assert o["n_neighbours"] == 9 and o["n_agreeing"] == 8
    # the sulfonamide cliff is real: CPD-013 and CPD-014 agree with each other, so neither is an outlier
    ids = {tuple(sorted(p["source_ids_a"] + p["source_ids_b"])): p["verdict"] for p in k["pairs"]}
    assert ids[("CPD-011", "CPD-013")] == "cliff" and ids[("CPD-012", "CPD-013")] == "cliff"


def test_split_and_generated(tmp_path):
    _run("split", EX / "bioactivity.csv", "-s", "scaffold", "-f", "train=0.8,test=0.2", "--seed", "0", "-o", tmp_path / "s.csv", "--json", tmp_path / "s.json")
    s = json.loads((tmp_path / "s.json").read_text())
    assert s["n_records"] == 33 and s["n_excluded"] == 3 and sum(s["n_by_split"].values()) == 30
    _run("split", EX / "bioactivity.csv", "-s", "temporal", "--date-column", "year", "-f", "train=0.7,test=0.3")
    _run("generated", EX / "generated.smi", "-r", EX / "train.smi", "--constraints", "mw=0:500,logp=-1:5", "--repair", "--json", tmp_path / "g.json")
    g = json.loads((tmp_path / "g.json").read_text())
    assert g["n_generated"] == 20 and g["n_valid"] == 18 and g["validity"] == pytest.approx(0.9)
    assert g["n_unique"] == 15 and g["novelty"] == pytest.approx(14 / 15, abs=1e-3)


def test_library_workflow(tmp_path):
    r = _run("validate", "--file", EX / "library.csv")
    assert "valid=21" in r.output and "invalid=3" in r.output
    _run("diagnose", "--file", EX / "library.csv", "-o", tmp_path / "d.csv")
    rows = list(csv.DictReader(open(tmp_path / "d.csv")))
    assert len(rows) == 3          # one row per failing record
    _run("standardize", "--file", EX / "library.csv", "-o", tmp_path / "s.csv")
    r = _run("identity", "--file", EX / "library.csv", "--level", "parent", "-o", tmp_path / "i.csv")
    assert "distinct at 'parent'=18" in r.output and "removable duplicates=3" in r.output
    r = _run("filter", "--file", EX / "library.csv", "--rules", "lipinski,pains", "--prep", "explicit-h", "-o", tmp_path / "f.csv")
    assert "Pass=18" in r.output and "Fail=3" in r.output and "Errors=3" in r.output
    r = _run("diff", EX / "library.csv", EX / "library_v2.csv", "--level", "parent", "--json", tmp_path / "diff.json")
    d = json.loads((tmp_path / "diff.json").read_text())
    assert (d["n_added"], d["n_removed"], d["n_unchanged"], d["n_changed"]) == (2, 2, 15, 1)
    assert d["changes_by_kind"] == {"salt, counter-ion or charge form": 1}
    for level, removed, changed in (("nostereo", 1, 2), ("skeleton", 0, 3)):
        _run("diff", EX / "library.csv", EX / "library_v2.csv", "--level", level, "--json", tmp_path / f"{level}.json")
        dd = json.loads((tmp_path / f"{level}.json").read_text())
        assert (dd["n_removed"], dd["n_changed"]) == (removed, changed), level


def test_quickstart_files(tmp_path):
    r = _run("validate", "--file", EX / "demo.csv")
    assert "valid=6" in r.output and "invalid=3" in r.output
    for level, added, removed, unchanged, changed in (("parent", 1, 2, 3, 0), ("nostereo", 1, 1, 2, 1), ("skeleton", 1, 0, 1, 2)):
        _run("diff", EX / "demo.csv", EX / "demo_v2.csv", "--level", level, "--json", tmp_path / f"{level}.json")
        d = json.loads((tmp_path / f"{level}.json").read_text())
        assert (d["n_added"], d["n_removed"], d["n_unchanged"], d["n_changed"]) == (added, removed, unchanged, changed), level


def test_search_and_structures(tmp_path):
    r = _run("similar", "CC(=O)Oc1ccccc1C(=O)O", "--file", EX / "library.csv", "--threshold", "0.5", "--top", "5")
    assert "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]" in r.output
    r = _run("substructure", "[NX3;H2]", "--file", EX / "library.csv")
    assert "Found 3 matches out of 24" in r.output
    r = _run("substructure", "c[H]", "--file", EX / "library.smi", "--prep", "explicit-h")
    assert "Found 14 matches out of 23" in r.output
    r = _run("rgroup", "c1ccccc1[*:1]", "--file", EX / "analogues.smi")
    assert r.output.count("Matched") == 6
    _run("scaffold", "--file", EX / "library.csv", "-o", tmp_path / "sc.csv")
    _run("fingerprint", "--file", EX / "library.csv", "--type", "ecfp4", "-o", tmp_path / "fp.csv")


def test_alert_catalogues_on_small_library(tmp_path):
    # library-independent defects, measured against the example library so the test stays fast
    _run("smartsaudit", EX / "alerts.csv", "--library", EX / "library.smi", "--checks", "all,proof", "--json", tmp_path / "a.json")
    a = json.loads((tmp_path / "a.json").read_text())
    st = {p["name"]: p for p in a["patterns"]}
    assert a["n_patterns"] == 18
    assert not st["unbalanced"]["parses"] and st["unbalanced"]["parse_error"]
    assert st["nitro_again"]["duplicate_of"] == 0 and st["nitro"]["index"] == 0
    assert st["aromatic_CH"]["requires_explicit_h"]
    assert st["oganesson"]["never_fires"] and st["oganesson"]["dead_verdict"] == "never-matching atom"
    assert st["recursive_carbonyl"]["has_recursive_smarts"] and st["recursive_carbonyl"]["proof_status"] == "not analysable"
    assert st["phenol"]["proven_subsumed_by"] is not None and st["benzene_kekule"]["proven_equivalent_to"] is not None
    _run("smartsaudit", EX / "alerts.smarts", "--library", EX / "library.smi")
    _run("smartsproof", EX / "alerts.csv", "--json", tmp_path / "p.json")
    p = json.loads((tmp_path / "p.json").read_text())
    assert p["n_patterns"] == 18 and p["n_proven_redundant"] == 12 and p["n_not_analysable"] >= 1
    _run("smartsdiff", EX / "alerts.csv", EX / "alerts_v2.csv", "--library", EX / "library.smi", "--json", tmp_path / "d.json")
    d = json.loads((tmp_path / "d.json").read_text())
    assert d["n_patterns_a"] == 18 and d["n_patterns_b"] == 13 and d["n_paired"] == 12
