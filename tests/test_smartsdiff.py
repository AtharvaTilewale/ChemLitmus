"""Tests for smartsdiff: pairing, semantic classification, verdict changes, CLI."""

import json

from typer.testing import CliRunner

from chemlitmus import diff_smarts
from chemlitmus.cli import app
from chemlitmus.core.smartsdiff import RDKIT_CATALOGS, load_side

runner = CliRunner()

LIB = ["CCO", "CC(=O)O", "c1ccccc1", "c1ccccc1O", "CCN", "C1CCCCC1", "OCCO", "CC(=O)Oc1ccccc1C(=O)O"]


def _write(tmp_path, name, rows):
    import csv
    p = tmp_path / name
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["name", "smarts"]); w.writerows(rows)
    return p


def _lib(tmp_path):
    p = tmp_path / "lib.smi"; p.write_text("\n".join(LIB) + "\n"); return p


def test_pairing_and_classification(tmp_path):
    a = _write(tmp_path, "a.csv", [("alcohol", "[OX2H]"), ("acid", "C(=O)[OH]"), ("amine", "[NX3;H2]"), ("ring6", "C1CCCCC1"),
                                   ("dead", "[Og]"), ("bad", "C(("), ("phenol_old", "c[OH]")])
    b = _write(tmp_path, "b.csv", [("alcohol", "[CX4][OX2H]"),          # narrowed: loses phenol
                                   ("acid", "C(=O)[OH]"),                # identical
                                   ("amine", "[NX3;H2,H1]"),             # text rewritten, same hits on this library
                                   ("ring6", "[R]"),                     # broadened: now also benzene/phenol
                                   ("dead", "[Og]"),
                                   ("bad", "C(=O)"),                     # repaired
                                   ("phenol_new", "c[OH]"),              # renamed -> paired by smarts
                                   ("ester", "C(=O)O[#6]")])             # added
    res = diff_smarts(a, b, library=_lib(tmp_path))
    by = {d.key: d for d in res.patterns}
    assert res.n_patterns_a == 7 and res.n_patterns_b == 8 and res.n_paired == 7
    assert res.paired_by == {"name": 6, "smarts": 1}
    assert by["alcohol"].semantic_status == "narrowed" and by["alcohol"].lost == 3        # acid OH, phenol, aspirin COOH
    assert set(by["alcohol"].examples_lost) == {"CC(=O)O", "Oc1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"}
    assert by["acid"].semantic_status == "same hits" and by["acid"].text_status == "identical"
    assert by["amine"].semantic_status == "same hits" and by["amine"].text_status == "rewritten"
    assert by["ring6"].semantic_status == "broadened" and by["ring6"].gained == 3 and by["ring6"].jaccard == 0.25
    assert by["dead"].semantic_status == "same hits" and by["dead"].hits_a == 0 and by["dead"].jaccard is None
    assert by["bad"].semantic_status == "repaired" and by["bad"].a.parses is False and by["bad"].b.parses
    assert by["phenol_old"].paired_by == "smarts" and by["phenol_old"].b.name == "phenol_new"
    assert by["ester"].semantic_status == "added" and by["ester"].hits_b == 1
    assert res.semantic_counts == {"same hits": 4, "narrowed": 1, "broadened": 1, "repaired": 1, "added": 1}
    # verdict change: A flags CCO, CC(=O)O, phenol, CCN, cyclohexane, OCCO, aspirin = 7; B flags all 8 (benzene via [R])
    assert res.flagged_a == 7 and res.flagged_b == 8 and res.flagged_only_b == 1 and res.flagged_only_a == 0
    assert res.verdict_changes == 1 and abs(res.verdict_change_fraction - 1 / 8) < 1e-9
    assert res.patterns[0].semantic_status == "repaired"              # changed patterns are listed first


def test_pair_by_identical_hits_and_removed(tmp_path):
    a = _write(tmp_path, "a.csv", [("x", "[OX2H][CX4]"), ("gone", "[Cl]")])
    b = _write(tmp_path, "b.csv", [("y", "[CX4][OX2H]")])                # renamed and rewritten: same hits
    res = diff_smarts(a, b, library=_lib(tmp_path))
    by = {d.key: d for d in res.patterns}
    assert by["x"].paired_by == "hits" and by["x"].semantic_status == "same hits" and by["x"].text_status == "rewritten"
    assert by["gone"].semantic_status == "removed"


def test_rdkit_catalog_side(tmp_path):
    sides, queries, cat = load_side("rdkit:PAINS_A")
    assert cat is not None and len(sides) == 16 and sides[0].smarts is None and queries == []
    assert "PAINS" in RDKIT_CATALOGS
    res = diff_smarts("rdkit:PAINS_A", "rdkit:PAINS_A", library=_lib(tmp_path))
    assert res.n_paired == 16 and res.semantic_counts == {"same hits": 16} and res.verdict_changes == 0
    assert res.text_counts == {"no SMARTS": 16}


def test_cli_smartsdiff(tmp_path):
    a = _write(tmp_path, "a.csv", [("alcohol", "[OX2H]")]); b = _write(tmp_path, "b.csv", [("alcohol", "[CX4][OX2H]")])
    out, js = tmp_path / "d.csv", tmp_path / "d.json"
    r = runner.invoke(app, ["smartsdiff", str(a), str(b), "--library", str(_lib(tmp_path)), "-o", str(out), "--json", str(js)])
    assert r.exit_code == 0, r.output
    res = diff_smarts(a, b, library=_lib(tmp_path))
    assert "narrowed" in r.output and f"Verdict changes: {res.verdict_changes} molecules" in r.output
    rows = out.read_text().splitlines(); assert rows[0].startswith("key,paired_by,semantic_status") and len(rows) == 2
    assert json.loads(js.read_text())["verdict_changes"] == res.verdict_changes == 3   # phenol, acetic acid, aspirin lose the flag
    assert runner.invoke(app, ["smartsdiff", str(a), str(b), "--prep", "nope"]).exit_code == 1
    assert runner.invoke(app, ["smartsdiff", str(a), "rdkit:NOPE"]).exit_code == 1
    assert runner.invoke(app, ["smartsdiff", str(a), str(tmp_path / "missing.csv")]).exit_code == 1
    r = runner.invoke(app, ["smartsdiff", "rdkit:PAINS_A", "rdkit:PAINS_A", "--library", str(_lib(tmp_path))])
    assert r.exit_code == 0 and "same hit set" in r.output
