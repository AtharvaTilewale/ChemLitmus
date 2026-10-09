"""Dataset audit, leakage, label conflicts, outputs and gates."""

import csv
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from chemlitmus import ChemicalPolicy, audit_dataset, label_conflicts, leakage_report
from chemlitmus.cli import app
from chemlitmus.core.dataset_audit import ISSUE_CATALOGUE
from chemlitmus.core.labels import parse_measurement
from chemlitmus.core.report import CleanPolicy, GatePolicy, evaluate_gates, write_audit_outputs

runner = CliRunner()

DS = """id,smiles,split,activity,units,relation,extra
A1,CC(=O)Oc1ccccc1C(=O)O,train,10,nM,=,x
A2,CC(=O)Oc1ccccc1C(=O)[O-].[Na+],test,20000,nM,=,y
A3,CC O,train,,,,
A4,,test,,,,
A5,c1ccccc1,train,1,uM,<,
A6,c1ccccc1,train,5,uM,=,
A7,C[C@H](N)C(=O)O,train,3,nM,=,
A8,CC(N)C(=O)O,test,4,nM,=,
A9,O=C1NC(=O)C=C1,test,,,,
A10,c1ccc2ccccc2c1,valid,7,,=,
A11,[13CH3]CO,train,1,mM,=,
A12,C(( ,test,,,,
A13,O=C(O)c1ccccc1O.O=C(O)c1ccccc1O,train,2,nM,=,
"""


@pytest.fixture
def ds(tmp_path):
    f = tmp_path / "ds.csv"; f.write_text(DS); return f


def _codes(audit, rid):
    return {i.code for i in audit.issues if i.record_id == rid}


def test_accounting_and_issue_codes(ds):
    a = audit_dataset(ds)
    s = a.summary
    assert (s.n_total, s.n_ok, s.n_empty, s.n_invalid, s.n_annotated) == (13, 10, 1, 2, 10)
    assert s.n_total == s.n_ok + s.n_empty + s.n_invalid + s.n_unsupported + s.n_error == len(a.records)
    assert s.processing_complete and s.identity_level == "parent"
    assert {"STD_CHANGED_IDENTITY", "DUP_PARENT", "SPLIT_OVERLAP", "LABEL_CONFLICT", "FRAG_MULTIPLE"} <= _codes(a, "r000001")   # aspirin sodium
    assert "DUP_EXACT" not in _codes(a, "r000001")
    assert _codes(a, "r000002") == {"PARSE_WHITESPACE", "REPAIR_CANDIDATE"} and _codes(a, "r000003") == {"PARSE_EMPTY"}
    assert "RELATION_CENSORED" in _codes(a, "r000004") and "DUP_EXACT" in _codes(a, "r000005")
    assert {"DUP_NOSTEREO", "SPLIT_RELATED", "STEREO_UNASSIGNED"} <= _codes(a, "r000007") and "SPLIT_OVERLAP" not in _codes(a, "r000007")
    assert "UNITS_MISSING" in _codes(a, "r000009") and "ISOTOPE" in _codes(a, "r000010") and "ELEMENT_UNUSUAL" not in _codes(a, "r000010")
    assert "PARSE_INVALID" in _codes(a, "r000011") and "FRAG_MULTIPLE_ORGANIC" in _codes(a, "r000012")
    for i in a.issues:
        assert i.code in ISSUE_CATALOGUE and i.severity in ("info", "warning", "error") and i.suggested_action
    # original fields survive
    assert a.records[0].fields["extra"] == "x" and a.records[0].source_id == "A1"
    # standardisation provenance
    tr = [t for t in a.records[1].transformations if t.changed_identity]
    assert tr and tr[0].operation == "largest_organic_fragment" and tr[0].identity_relation.startswith("salt")
    assert a.records[1].standardized_smiles == a.records[0].standardized_smiles
    # identity groups
    g = next(g for g in a.groups if g.size == 2 and "r000001" in g.record_ids)
    assert g.level == "parent" and g.splits == ["test", "train"] and "salt, counter-ion or charge form" in g.differs_by
    assert s.n_groups_by_level["nostereo"] < s.n_groups_by_level["parent"]


def test_leakage_classes_are_distinct_and_formula_is_not_identity(ds):
    a = audit_dataset(ds)
    pair = next(p for p in a.leakage["pairs"] if p["evaluation_split"] == "test")
    ov = {o["level"]: o["n_eval_records"] for o in pair["overlap"]}
    assert ov == {"exact": 0, "parent": 1, "tautomer": 1, "nostereo": 2, "skeleton": 2}
    assert pair["formula_matches"] == 2 and pair["n_related_by_similarity"] >= 1 and pair["neighbour_fingerprint"] == "ecfp4/2048"
    assert a.leakage["within_split_duplicates"]["train"] == 1                      # benzene twice
    assert "not additive" in a.leakage["note"]
    # standalone API with curated salt / tautomer / stereoisomer / constitutional-isomer cases
    rep = leakage_report({
        "train": [("t1", "CC(=O)Oc1ccccc1C(=O)O", None), ("t2", "Oc1ccccn1", None), ("t3", "C[C@H](N)C(=O)O", None), ("t4", "CCCO", None)],
        "test": [("e1", "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]", None), ("e2", "O=c1cccc[nH]1", None), ("e3", "C[C@@H](N)C(=O)O", None), ("e4", "CC(C)O", None)],
    })
    ov = {o.level: {ex["eval_record_id"] for ex in o.examples} for o in rep.pairs[0].overlap}
    assert ov["exact"] == set() and ov["parent"] == {"e1"} and "e2" in ov["tautomer"] and "e3" in ov["nostereo"] and {"e1", "e2", "e3"} <= ov["skeleton"]
    assert "e4" not in ov["skeleton"] and rep.pairs[0].formula_matches >= 1          # propanol isomers: formula only
    assert rep.pairs[0].overlap[1].n_eval_records == 1


def test_temporal_and_reference_selection():
    rep = leakage_report({"a": [("1", "CCO", "2020-01-01"), ("2", "CCN", "2021-06-01")], "b": [("3", "CCCO", "2021-01-01"), ("4", "CCCN", "2022-01-01")]}, reference="a")
    p = rep.pairs[0]
    assert p.reference_split == "a" and p.temporal_violations == 1 and "3" in rep.record_issues
    rep2 = leakage_report({"small": [("1", "CCO", None)], "big": [("2", "CCN", None), ("3", "CCC", None)]})
    assert rep2.pairs[0].reference_split == "big"


def test_label_conflicts_units_relations_context():
    recs = [
        {"record_id": "1", "key": "K", "val": "10", "u": "nM", "rel": "=", "target": "T1"},
        {"record_id": "2", "key": "K", "val": "20", "u": "uM", "rel": "=", "target": "T1"},     # 2000-fold -> conflict (3.3 log)
        {"record_id": "3", "key": "K", "val": "1", "u": "uM", "rel": ">", "target": "T1"},      # censored: excluded from spread
        {"record_id": "4", "key": "K", "val": "15", "u": "nM", "rel": "=", "target": "T2"},     # different target: separate context
        {"record_id": "5", "key": "K", "val": "12", "u": "", "rel": "=", "target": "T2"},       # missing units
        {"record_id": "6", "key": "L", "val": "5", "u": "nM", "rel": "=", "target": "T1"},
        {"record_id": "7", "key": "L", "val": "6", "u": "nM", "rel": "=", "target": "T1"},      # within tolerance
    ]
    rep = label_conflicts(recs, endpoint_field="val", units_field="u", relation_field="rel", context_fields=["target"])
    assert rep.kind == "quantitative" and rep.n_groups_compared == 3 and rep.n_conflicts == 1
    g = rep.groups[0]
    assert g.context == {"target": "T1"} and g.n_exact == 2 and g.n_censored == 1 and g.scale == "log10 molar" and abs(g.spread - 3.301) < 0.01
    assert rep.n_censored == 1 and rep.n_missing_units == 1 and "5" in rep.record_issues and rep.record_issues["5"][0]["code"] == "UNITS_MISSING"
    assert all(d["code"] != "LABEL_CONFLICT" for d in rep.record_issues.get("4", []))
    m = parse_measurement("<5", "uM")
    assert m.censored and m.relation == "<" and m.value_nm == 5000 and abs(m.log_value - 5.301) < 0.01 and "-> 5000 nM" in m.conversion
    assert parse_measurement("1-10", "nM").relation == "range" and parse_measurement("abc", "nM").comparable is False
    rep_c = label_conflicts([{"record_id": "1", "key": "K", "y": "active"}, {"record_id": "2", "key": "K", "y": "inactive"}, {"record_id": "3", "key": "L", "y": "active"}], endpoint_field="y")
    assert rep_c.kind == "classification" and rep_c.n_conflicts == 1 and rep_c.groups[0].labels == {"active": 1, "inactive": 1}
    # salts group together at parent level only when the keys say so: the caller supplies keys, so check mixed-unit incomparability
    rep_m = label_conflicts([{"record_id": "1", "key": "K", "val": "10", "u": "nM"}, {"record_id": "2", "key": "K", "val": "10", "u": "percent"}], endpoint_field="val", units_field="u")
    assert rep_m.n_conflicts == 0 and rep_m.n_groups_compared == 1


def test_outputs_gates_and_clean_export(ds, tmp_path):
    a = audit_dataset(ds)
    out = tmp_path / "out"
    paths = write_audit_outputs(a, out, input_path=ds, gates=GatePolicy(no_split_overlap=True), clean=CleanPolicy(), command="test")
    assert set(Path(p).name for p in paths.values()) >= {"records.csv", "issues.csv", "identity_groups.csv", "summary.json", "audit.json", "policy.json", "report.html", "manifest.json", "clean.csv", "exclusions.csv"}
    rows = list(csv.DictReader(open(paths["records"])))
    assert len(rows) == 13 and rows[0]["in:extra"] == "x" and rows[1]["identity_parent"] == rows[0]["identity_parent"] and rows[2]["repair_candidate"] == "CCO"
    summary = json.loads(Path(paths["summary"]).read_text())
    assert summary["gate"]["passed"] is False and summary["gate"]["exit_code"] == 3 and summary["policy_hash"] == a.policy_hash
    clean = list(csv.DictReader(open(paths["clean"]))); excl = list(csv.DictReader(open(paths["exclusions"])))
    assert len(clean) + len(excl) == 13
    assert {e["record_id"] for e in excl} == {"r000002", "r000003", "r000005", "r000011"}                 # whitespace, empty, exact dup, invalid
    assert all(c["structure_source"] == "standardized" for c in clean) and clean[0]["original_structure"] == "CC(=O)Oc1ccccc1C(=O)O"
    man = json.loads(Path(paths["manifest"]).read_text())
    assert man["inputs"][0]["path"] == "ds.csv" and man["policy_hash"] == a.policy_hash and len(man["outputs"]) == 9 and man["rdkit_version"]
    html = Path(paths["report"]).read_text()
    assert "<svg" in html and "SPLIT_OVERLAP" in html and "not additive" in html and "&lt;" not in html[:200]
    # gates
    assert evaluate_gates(a, GatePolicy()).exit_code == 0
    assert evaluate_gates(a, GatePolicy(fail_on_severity="error")).exit_code == 3
    assert evaluate_gates(a, GatePolicy(max_invalid_fraction=0.1)).exit_code == 3 and evaluate_gates(a, GatePolicy(max_invalid_fraction=0.5)).exit_code == 0
    assert evaluate_gates(a, GatePolicy(no_label_conflicts=True)).exit_code == 3
    assert evaluate_gates(a, GatePolicy(max_duplicate_fraction=0.0)).exit_code == 3


def test_policy_changes_behaviour(ds):
    cons = audit_dataset(ds, policy="conservative")
    assert cons.summary.identity_level == "exact" and not any(i.code == "STD_CHANGED_IDENTITY" for i in cons.issues)
    assert not any(i.code == "REPAIR_CANDIDATE" for i in cons.issues)              # repair.mode = none
    assert not any(i.code == "SPLIT_OVERLAP" for i in cons.issues)                 # salt vs acid are distinct at 'exact'
    pol = ChemicalPolicy(identity_level="nostereo", severity={"overrides": {"ALERT_MATCH": "ignore", "DUP_NOSTEREO": "error"}})
    a = audit_dataset(ds, policy=pol)
    assert not any(i.code == "ALERT_MATCH" for i in a.issues)
    assert any(i.code == "SPLIT_OVERLAP" and i.record_id == "r000007" for i in a.issues)   # alanine stereo forms overlap at nostereo


def test_html_escapes_user_text(tmp_path):
    f = tmp_path / "x.csv"; f.write_text('id,smiles\n"<script>alert(1)</script>",CC O\nB,C((\n')
    a = audit_dataset(f)
    from chemlitmus.core.report import render_html
    html = render_html(a)
    assert "<script>alert(1)</script>" not in html and "&lt;script&gt;" in html


def test_cli_audit_leakage_conflicts(ds, tmp_path):
    out = tmp_path / "o"
    r = runner.invoke(app, ["audit", str(ds), "-o", str(out), "--clean", "--no-split-overlap"])
    assert r.exit_code == 3, r.output
    assert "Gate: FAIL" in r.output and (out / "report.html").exists() and (out / "clean.csv").exists()
    r = runner.invoke(app, ["audit", str(ds), "-o", str(out), "--config", "conservative"])
    assert r.exit_code == 0 and "Gate: PASS" in r.output
    assert runner.invoke(app, ["audit", str(tmp_path / "missing.csv")]).exit_code == 1
    assert runner.invoke(app, ["audit", str(ds), "--fail-on-severity", "high"]).exit_code == 1
    amb = tmp_path / "amb.csv"; amb.write_text("foo,bar\nCCO,1\n")
    assert runner.invoke(app, ["audit", str(amb), "-o", str(out)]).exit_code == 1
    assert runner.invoke(app, ["audit", str(amb), "-o", str(out), "--first-column"]).exit_code == 0
    r = runner.invoke(app, ["leakage", str(ds), "--split-column", "split", "--fail-on-overlap", "--json", str(tmp_path / "l.json"), "-o", str(tmp_path / "l.csv")])
    assert r.exit_code == 3 and "policy level" in r.output and json.loads((tmp_path / "l.json").read_text())["identity_level"] == "parent"
    tr = tmp_path / "train.csv"; te = tmp_path / "test.csv"
    tr.write_text("smiles\nCCO\nc1ccccc1\n"); te.write_text("smiles\nCCO\nCCN\n")
    r = runner.invoke(app, ["leakage", str(tr), str(te)])
    assert r.exit_code == 0 and "train → test" in r.output
    assert runner.invoke(app, ["leakage", str(tr)]).exit_code == 1
    r = runner.invoke(app, ["conflicts", str(ds), "-e", "activity", "--units-column", "units", "--relation-column", "relation", "--json", str(tmp_path / "c.json")])
    assert r.exit_code == 0 and "1 in conflict" in r.output
    assert runner.invoke(app, ["conflicts", str(ds), "-e", "activity", "--kind", "weird"]).exit_code == 1
    assert runner.invoke(app, ["conflicts", str(ds), "-e", "activity", "--context", "nope"]).exit_code == 1
