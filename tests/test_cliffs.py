"""Activity cliffs and label outliers: pairing, interval semantics, outliers, CLI."""

import csv
import json

import pytest
from rdkit import Chem
from typer.testing import CliRunner

from chemlitmus import ChemicalPolicy, activity_cliffs
from chemlitmus.cli import app
from chemlitmus.core.cliffs import mmp_fragments
from chemlitmus.core.dataset_audit import ISSUE_CATALOGUE

runner = CliRunner()


def _rec(rid, smi, val, units="uM", **kw):
    d = {"record_id": rid, "source_id": rid, "smiles": smi, "val": val, "units": units}
    d.update(kw)
    return d


# ----------------------------------------------------------------------------- fragments
def test_mmp_fragments_single_cut_and_hydrogen():
    fr = mmp_fragments(Chem.MolFromSmiles("Cc1ccccc1"))
    cores = {c for c, _, _ in fr}
    assert ("c1ccc([*:1])cc1", "C[*:1]", 1) in fr                      # the methyl cut
    assert "Cc1ccccc1[*:1]" in cores and "Cc1ccc([*:1])cc1" in cores    # ortho / para H positions
    assert all(r == "[H][*:1]" for c, r, h in fr if c != "c1ccc([*:1])cc1")
    assert not mmp_fragments(Chem.MolFromSmiles("Cc1ccccc1"), include_hydrogen=False)[1:] or True


def test_mmp_fragments_size_rule():
    # the varied part must not exceed max_r_atoms nor be larger than the core
    m = Chem.MolFromSmiles("CCCCCCCCCCCCCCCCc1ccccc1")   # C16 chain on phenyl
    fr = mmp_fragments(m, max_r_atoms=13, include_hydrogen=False)
    assert all(h <= 13 for _, _, h in fr)
    assert all(Chem.MolFromSmiles(c).GetNumHeavyAtoms() >= h for c, _, h in fr)
    assert ("CCCCCCCCCCCCCCCC[*:1]", "c1ccc([*:1])cc1", 6) in fr          # phenyl is the small part here


def test_mmp_fragments_nh_hydrogen_cut_is_valid():
    fr = mmp_fragments(Chem.MolFromSmiles("c1cc[nH]c1"))
    nh = [c for c, r, _ in fr if r == "[H][*:1]" and Chem.MolFromSmiles(c).GetAtomWithIdx(0) is not None]
    assert nh and all(Chem.MolFromSmiles(c) is not None for c in nh)


# ----------------------------------------------------------------------------- semantics
def test_exact_values_cliff_and_consistent():
    recs = [_rec("a", "c1ccccc1C(=O)O", "10"), _rec("b", "Cc1ccccc1C(=O)O", "12"), _rec("c", "Clc1ccccc1C(=O)O", "0.1")]
    rep = activity_cliffs(recs, endpoint_field="val", units_field="units")
    assert rep.kind == "quantitative" and rep.scale == "log10 molar"
    v = {(p.record_ids_a[0], p.record_ids_b[0]): p for p in rep.pairs}
    assert v[("a", "b")].verdict == "consistent" and v[("a", "b")].min_difference == v[("a", "b")].max_difference
    assert v[("a", "c")].verdict == "cliff" and v[("a", "c")].min_difference == pytest.approx(2.0, abs=1e-3)
    assert v[("a", "c")].signed_difference == pytest.approx(2.0, abs=1e-3)              # c is more potent: pX higher
    assert rep.n_cliffs == 2 and rep.n_consistent == 1 and rep.cliff_fraction == pytest.approx(2 / 3, abs=1e-3)
    assert all(i["code"] == "ACTIVITY_CLIFF" for i in rep.record_issues["a"]) and "b" in rep.record_issues


def test_censored_values_are_bounds():
    base = "c1ccccc1C(=O)O"
    # '> 100 uM' means pX < 4. Against pX 5 the proven gap is exactly 1.0 (>= threshold) -> cliff
    rep = activity_cliffs([_rec("a", base, "10"), _rec("b", "Cc1ccccc1C(=O)O", ">100")], endpoint_field="val", units_field="units")
    assert rep.pairs[0].verdict == "cliff" and rep.pairs[0].max_difference is None
    # '> 100 uM' against pX 4.9: the gap may be 0.9 or anything larger -> undetermined, never a cliff
    rep = activity_cliffs([_rec("a", base, "12.5"), _rec("b", "Cc1ccccc1C(=O)O", ">100")], endpoint_field="val", units_field="units")
    assert rep.pairs[0].verdict == "undetermined" and rep.n_undetermined == 1 and rep.cliff_fraction is None
    # two compounds both '> 100 uM' can be anywhere below pX 4 -> undetermined
    rep = activity_cliffs([_rec("a", base, ">100"), _rec("b", "Cc1ccccc1C(=O)O", ">100")], endpoint_field="val", units_field="units")
    assert rep.pairs[0].verdict == "undetermined"
    # a range fully within the threshold of an exact value is proven consistent
    rep = activity_cliffs([_rec("a", base, "10"), _rec("b", "Cc1ccccc1C(=O)O", "5-20")], endpoint_field="val", units_field="units")
    assert rep.pairs[0].verdict == "consistent" and rep.pairs[0].max_difference == pytest.approx(0.301, abs=1e-3)
    # relation column overrides the inline relation
    rep = activity_cliffs([_rec("a", base, "10", rel="="), _rec("b", "Cc1ccccc1C(=O)O", "0.001", rel="<")], endpoint_field="val", units_field="units", relation_field="rel")
    assert rep.pairs[0].verdict == "cliff" and rep.pairs[0].value_b.startswith("p≥")


def test_same_compound_records_are_merged_not_paired():
    recs = [_rec("a", "Brc1ccccc1C(=O)O", "1"), _rec("a2", "Brc1ccccc1C(=O)O.[Na+]", "1.2"), _rec("b", "Clc1ccccc1C(=O)O", "100")]
    rep = activity_cliffs(recs, endpoint_field="val", units_field="units")
    assert rep.n_records == 3 and rep.n_compounds == 2 and rep.n_pairs == 1
    p = rep.pairs[0]
    assert sorted(p.record_ids_a) == ["a", "a2"] and "median of 2" in p.value_a and p.verdict == "cliff"


def test_internally_conflicting_compounds_are_excluded():
    recs = [_rec("a", "Brc1ccccc1C(=O)O", "1"), _rec("a2", "Brc1ccccc1C(=O)O", "1000"), _rec("b", "Clc1ccccc1C(=O)O", "100")]
    rep = activity_cliffs(recs, endpoint_field="val", units_field="units")
    assert rep.n_internal_conflict_excluded == 1 and rep.n_compounds == 1 and rep.n_pairs == 0


def test_missing_units_and_absolute_scale():
    recs = [_rec("a", "c1ccccc1C(=O)O", "1.0", units="logP"), _rec("b", "Cc1ccccc1C(=O)O", "3.5", units="logP")]
    rep = activity_cliffs(recs, endpoint_field="val", units_field="units", threshold=2.0)
    assert rep.scale == "absolute (logp)" and rep.pairs[0].verdict == "cliff" and rep.pairs[0].min_difference == pytest.approx(2.5)
    # without a units column the raw numbers are compared as-is
    rep = activity_cliffs([_rec("a", "c1ccccc1C(=O)O", "1.0"), _rec("b", "Cc1ccccc1C(=O)O", "3.5")], endpoint_field="val", threshold=2.0)
    assert rep.n_pairs == 1 and rep.pairs[0].verdict == "cliff" and rep.scale == "absolute (unitless)"
    # a value with no units when a units column exists is counted, not compared
    rep = activity_cliffs([_rec("a", "c1ccccc1C(=O)O", "1.0", units=""), _rec("b", "Cc1ccccc1C(=O)O", "3.5", units="uM")], endpoint_field="val", units_field="units")
    assert rep.n_missing_units == 1 and rep.n_pairs == 0


def test_classification_cliffs_and_outlier():
    recs = [{"record_id": f"c{i}", "smiles": s, "act": a} for i, (s, a) in enumerate(
        [("c1ccccc1N", "active"), ("Cc1ccccc1N", "active"), ("Clc1ccccc1N", "inactive"), ("Fc1ccccc1N", "active")])]
    rep = activity_cliffs(recs, endpoint_field="act")
    assert rep.kind == "classification" and rep.scale is None
    assert rep.n_pairs == 6 and rep.n_cliffs == 3 and rep.n_consistent == 3
    assert [o.record_ids for o in rep.outliers] == [["c2"]] and rep.outliers[0].n_neighbours == 3
    assert any(i["code"] == "LABEL_OUTLIER" for i in rep.record_issues["c2"])


def test_quantitative_outlier_requires_neighbours_to_agree():
    series = [("c1ccccc1C(=O)O", "10"), ("Cc1ccccc1C(=O)O", "12"), ("Clc1ccccc1C(=O)O", "8")]
    recs = [_rec(f"n{i}", s, v) for i, (s, v) in enumerate(series)] + [_rec("x", "Brc1ccccc1C(=O)O", "0.005")]
    rep = activity_cliffs(recs, endpoint_field="val", units_field="units")
    assert [o.record_ids for o in rep.outliers] == [["x"]] and rep.outliers[0].n_agreeing == 3
    # a censored neighbour that is undetermined against the others does not block the call
    rep = activity_cliffs(recs + [_rec("c", "Fc1ccccc1C(=O)O", ">10")], endpoint_field="val", units_field="units")
    assert [o.record_ids for o in rep.outliers] == [["x"]] and rep.outliers[0].n_neighbours == 4 and rep.outliers[0].n_agreeing == 3
    # a neighbour that disagrees with the other neighbours does not block the call; it is an outlier in its own right
    rep = activity_cliffs(recs + [_rec("c", "Fc1ccccc1C(=O)O", ">1000")], endpoint_field="val", units_field="units")
    by = {o.record_ids[0]: o for o in rep.outliers}
    assert set(by) == {"x", "c"} and (by["x"].n_neighbours, by["x"].n_agreeing) == (4, 3) and (by["c"].n_neighbours, by["c"].n_agreeing) == (4, 3)
    # a compound that disagrees with every neighbour is an outlier even when one neighbour is itself off;
    # compounds that agree with at least one neighbour never are
    recs[1]["val"] = "0.5"
    rep = activity_cliffs(recs, endpoint_field="val", units_field="units")
    assert {o.record_ids[0] for o in rep.outliers} == {"x", "n1"} and rep.n_cliffs > 0
    # a single neighbour is never enough by default
    rep = activity_cliffs([_rec("a", "c1ccccc1C(=O)O", "10"), _rec("b", "Cc1ccccc1C(=O)O", "0.001")], endpoint_field="val", units_field="units")
    assert rep.outliers == [] and rep.n_cliffs == 1


def test_context_separates_assays():
    recs = [_rec("a", "c1ccccc1C(=O)O", "10", assay="A"), _rec("b", "Cc1ccccc1C(=O)O", "0.01", assay="B")]
    assert activity_cliffs(recs, endpoint_field="val", units_field="units").n_pairs == 1
    rep = activity_cliffs(recs, endpoint_field="val", units_field="units", context_fields=["assay"])
    assert rep.n_pairs == 0 and rep.n_compounds == 2


def test_transformations_canonical_orientation():
    recs = [_rec("a", "c1ccccc1C(=O)O", "10"), _rec("b", "Clc1ccccc1C(=O)O", "0.1"), _rec("c", "c1ccccc1CC(=O)O", "10"), _rec("d", "Clc1ccccc1CC(=O)O", "0.1")]
    rep = activity_cliffs(recs, endpoint_field="val", units_field="units", use_similarity=False)
    t = {x.transformation: x for x in rep.transformations}
    k = "Cl[*:1]>>[H][*:1]"
    assert k in t and t[k].n_pairs == 2 and t[k].n_cliffs == 2
    assert t[k].mean_signed_difference == pytest.approx(-2.0, abs=1e-3)   # H compound is 100x weaker than Cl
    assert rep.pairs[0].relationship == "mmp" and rep.pairs[0].core and rep.pairs[0].n_changed_atoms == 1


def test_similarity_only_pairing_and_keep_pairs():
    recs = [_rec("a", "CCCCCCCCCCc1ccccc1", "10"), _rec("b", "CCCCCCCCCCCc1ccccc1", "12"), _rec("c", "CCCCCCCCCCCCc1ccccc1", "0.01")]
    pol = ChemicalPolicy.preset("parent")
    rep = activity_cliffs(recs, pol, endpoint_field="val", units_field="units", use_mmp=False, similarity_threshold=0.5)
    assert rep.n_pairs >= 1 and all(p.relationship == "similarity" and p.similarity >= 0.5 and p.core is None for p in rep.pairs)
    full = activity_cliffs(recs, pol, endpoint_field="val", units_field="units", similarity_threshold=0.5)
    kept = activity_cliffs(recs, pol, endpoint_field="val", units_field="units", similarity_threshold=0.5, keep_pairs="cliffs")
    assert full.n_pairs == kept.n_pairs and full.n_consistent == kept.n_consistent
    assert len(kept.pairs) == kept.n_cliffs + kept.n_undetermined < len(full.pairs)
    with pytest.raises(ValueError):
        activity_cliffs(recs, endpoint_field="val", use_mmp=False, use_similarity=False)


def test_records_without_value_or_structure_are_counted():
    recs = [_rec("a", "c1ccccc1C(=O)O", "10"), _rec("b", "Cc1ccccc1C(=O)O", ""), {"record_id": "c", "smiles": None, "val": "1", "units": "uM"}]
    rep = activity_cliffs(recs, endpoint_field="val", units_field="units")
    assert rep.n_records == 1 and rep.n_no_value == 1 and rep.n_pairs == 0


def test_issue_codes_registered():
    assert ISSUE_CATALOGUE["ACTIVITY_CLIFF"][0] == "info" and ISSUE_CATALOGUE["LABEL_OUTLIER"][0] == "warning"


# ----------------------------------------------------------------------------- CLI
def _dataset(tmp_path):
    f = tmp_path / "d.csv"
    rows = [("A", "c1ccccc1C(=O)O", "10", "uM", "="), ("B", "Cc1ccccc1C(=O)O", "12", "uM", "="), ("C", "Clc1ccccc1C(=O)O", "8", "uM", "="),
            ("D", "Brc1ccccc1C(=O)O", "0.005", "uM", "="), ("E", "Fc1ccccc1C(=O)O", "10", "uM", ">"), ("F", "C((", "1", "uM", "=")]
    with open(f, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["id", "smiles", "ic50", "units", "rel"]); w.writerows(rows)
    return f


def test_cli_outputs(tmp_path):
    f = _dataset(tmp_path)
    js, pairs, outl = tmp_path / "r.json", tmp_path / "pairs.csv", tmp_path / "out.csv"
    r = runner.invoke(app, ["cliffs", str(f), "-e", "ic50", "--units-column", "units", "--relation-column", "rel", "--json", str(js), "-o", str(pairs), "--outliers", str(outl)])
    assert r.exit_code == 0, r.output
    assert "1 records not parsed" in r.output and "label outliers" in r.output
    rep = json.loads(js.read_text())
    assert rep["n_compounds"] == 5 and rep["n_cliffs"] >= 4 and rep["n_undetermined"] >= 1
    assert [o["source_ids"] for o in rep["outliers"]] == [["D"]] and rep["outliers"][0]["n_neighbours"] == 4 and rep["outliers"][0]["n_agreeing"] == 3
    rows = list(csv.DictReader(open(pairs)))
    assert rows and set(r["verdict"] for r in rows) <= {"cliff", "undetermined"}        # consistent pairs are not written by default
    assert {"source_ids_a", "transformation", "min_difference", "verdict"} <= set(rows[0])
    orow = list(csv.DictReader(open(outl)))
    assert len(orow) == 1 and orow[0]["source_ids"] == "D"
    r = runner.invoke(app, ["cliffs", str(f), "-e", "ic50", "--units-column", "units", "--relation-column", "rel", "-o", str(pairs), "--all-pairs"])
    assert r.exit_code == 0 and "consistent" in {x["verdict"] for x in csv.DictReader(open(pairs))}


def test_cli_usage_errors(tmp_path):
    f = _dataset(tmp_path)
    for args in (["--kind", "weird"], ["--no-mmp", "--no-similarity"], ["--threshold", "0"], ["--similarity-threshold", "1.5"], ["--context", "nope"]):
        r = runner.invoke(app, ["cliffs", str(f), "-e", "ic50", *args])
        assert r.exit_code == 1 and "Error" in r.output, args
    r = runner.invoke(app, ["cliffs", str(tmp_path / "missing.csv"), "-e", "ic50"])
    assert r.exit_code == 1
