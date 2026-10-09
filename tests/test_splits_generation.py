"""Group-aware split generation and generated-molecule evaluation."""

import csv
import json

import pytest
from typer.testing import CliRunner

from chemlitmus import evaluate_generated, make_splits
from chemlitmus.cli import app

runner = CliRunner()

SMILES = ["CCO", "CCO.[Na]", "c1ccccc1", "c1ccccc1C", "c1ccccc1CC", "CCN", "CCCN", "c1ccncc1", "CCCCO", "C[C@H](N)C(=O)O"]
RECS = [{"record_id": f"r{i}", "source_id": f"S{i}", "smiles": s, "y": "a" if i % 2 else "b",
         "date": f"20{10 + i}-01-01", "source": f"doc{i // 3}"} for i, s in enumerate(SMILES)]


def test_identity_groups_are_indivisible_and_seeded():
    r = make_splits(RECS, {"train": 0.7, "test": 0.3}, strategy="identity", seed=1, endpoint_field="y")
    by_rid = {a.record_id: a for a in r.assignments}
    assert by_rid["r0"].split == by_rid["r1"].split                          # ethanol and its sodium salt share a parent key
    assert by_rid["r0"].group_key == by_rid["r1"].group_key
    assert sum(r.n_by_split.values()) == len(RECS) and r.n_groups == 9
    assert make_splits(RECS, {"train": 0.7, "test": 0.3}, strategy="identity", seed=1).n_by_split == r.n_by_split
    other = make_splits(RECS, {"train": 0.7, "test": 0.3}, strategy="identity", seed=7)
    assert {a.record_id: a.split for a in other.assignments} != {a.record_id: a.split for a in r.assignments} or other.n_by_split == r.n_by_split
    assert set(r.endpoint_balance) == {"train", "test"}
    ov = {o["level"]: o["n_eval_records"] for o in r.leakage["pairs"][0]["overlap"]}
    assert ov["parent"] == 0 and ov["exact"] == 0                            # verified with the audit's own implementation


def test_scaffold_split_does_not_pool_acyclic_molecules():
    r = make_splits(RECS, {"train": 0.6, "test": 0.4}, strategy="scaffold", seed=0)
    groups = {a.record_id: a.group_key for a in r.assignments}
    acyclic = [g for rid, g in groups.items() if g.startswith("acyclic:")]
    assert len(acyclic) == len(set(acyclic)) == 6                            # every acyclic record is its own group
    benzenes = {groups["r2"], groups["r3"], groups["r4"]}
    assert len(benzenes) == 1                                                # same Murcko scaffold -> one group
    assert r.n_groups == 8 and sum(r.n_by_split.values()) == len(RECS)


def test_impossible_constraints_are_reported_not_hidden():
    big = [{"record_id": f"b{i}", "smiles": "CCO"} for i in range(8)] + [{"record_id": "x", "smiles": "c1ccccc1"}]
    r = make_splits(big, {"train": 0.5, "test": 0.5}, strategy="identity", seed=0)
    assert r.conflicts and any("largest group" in c for c in r.conflicts) and any("achieved" in c for c in r.conflicts)
    assert max(r.n_by_split.values()) == 8                                   # the group was not broken
    with pytest.raises(ValueError):
        make_splits(RECS, {"train": 0.5, "test": 0.4})
    with pytest.raises(ValueError):
        make_splits(RECS, strategy="nope")
    bad = [{"record_id": "a", "smiles": "C(("}, {"record_id": "b", "smiles": "CCO"}]
    r2 = make_splits(bad, {"train": 0.5, "test": 0.5}, strategy="identity")
    assert r2.n_excluded == 1 and any("excluded" in c for c in r2.conflicts)


def test_temporal_and_source_strategies():
    r = make_splits(RECS, {"train": 0.7, "test": 0.3}, strategy="temporal", endpoint_field="y")
    dates = {a.split: sorted(next(x["date"] for x in RECS if x["record_id"] == a.record_id) for a in r.assignments if a.split == s) for s in ("train", "test") for a in r.assignments}
    train_max = max(x["date"] for x in RECS if any(a.record_id == x["record_id"] and a.split == "train" for a in r.assignments))
    test_min = min(x["date"] for x in RECS if any(a.record_id == x["record_id"] and a.split == "test" for a in r.assignments))
    assert train_max < test_min and r.seed is None
    del dates
    rs = make_splits(RECS, {"train": 0.7, "test": 0.3}, strategy="source")
    groups = {a.group_key for a in rs.assignments if a.split == "test"}
    assert groups and not (groups & {a.group_key for a in rs.assignments if a.split == "train"})   # a source is never split
    missing = [{"record_id": "a", "smiles": "CCO"}]
    assert make_splits(missing, strategy="temporal").n_excluded == 1


def test_generation_metrics_and_denominators():
    g = evaluate_generated(["CCO", "CCO", "OCC", "c1ccccc1", "bad((", "", "CC O", "CC(=O)Oc1ccccc1C(=O)O"],
                           reference={"train": ["CCO", "CCCC"]}, constraints={"mw": (0, 200)}, repair=True)
    assert g.n_generated == 8 and g.n_valid == 5 and g.validity == 0.625
    assert g.invalid_reasons == {"unparseable": 1, "empty": 1, "whitespace": 1} and g.n_empty == 1
    assert g.n_unique == 3 and g.uniqueness == 0.6                           # CCO thrice (CCO, CCO, OCC)
    assert g.n_novel == 2 and g.novelty == pytest.approx(2 / 3, abs=1e-3) and g.reference_sets == ["train"] and g.n_reference == 2
    assert g.molecules[0].novel is False and g.molecules[2].duplicate_of == 0
    assert g.nearest_neighbour_similarity["max"] == 1.0 and g.fingerprint == "ecfp4/2048"
    assert g.constraints == {"mw": 5} and g.scaffold_diversity is not None
    assert g.n_repaired == 1 and g.repaired_report is not None and g.repaired_report.n_valid == 1
    assert "aspirin" not in str(g.alerts_by_set)                             # alerts are descriptions, not verdicts
    assert all(lvl in g.uniqueness_by_level for lvl in ("exact", "parent", "skeleton"))


def test_generation_salt_and_stereo_novelty_follow_policy():
    from chemlitmus import ChemicalPolicy
    gen = ["CC(=O)Oc1ccccc1C(=O)[O-].[Na+]", "C[C@@H](N)C(=O)O"]
    ref = {"train": ["CC(=O)Oc1ccccc1C(=O)O", "C[C@H](N)C(=O)O"]}
    exact = evaluate_generated(gen, reference=ref, policy=ChemicalPolicy(identity_level="exact"))
    parent = evaluate_generated(gen, reference=ref, policy=ChemicalPolicy(identity_level="parent"))
    nostereo = evaluate_generated(gen, reference=ref, policy=ChemicalPolicy(identity_level="nostereo"))
    assert exact.novelty == 1.0                                              # both differ exactly
    assert parent.novelty == 0.5 and parent.molecules[0].novel is False      # the salt is not a new compound
    assert nostereo.novelty == 0.0                                           # neither salt nor enantiomer is new


def test_empty_valid_population_returns_undefined_not_zero():
    g = evaluate_generated(["bad((", "", "   "])
    assert g.n_valid == 0 and g.validity == 0.0 and g.novelty is None
    assert set(g.undefined) == {"uniqueness", "novelty", "scaffold_diversity", "nearest_neighbour_similarity"}
    assert g.descriptor_summary == {} and g.scaffold_diversity is None
    assert evaluate_generated([]).n_generated == 0


def test_cli_split_and_generated(tmp_path):
    f = tmp_path / "d.csv"
    f.write_text("id,smiles,y,date,source\n" + "\n".join(f"S{i},{s},{'a' if i % 2 else 'b'},20{10 + i}-01-01,doc{i // 3}" for i, s in enumerate(SMILES)) + "\n")
    out, js = tmp_path / "s.csv", tmp_path / "s.json"
    r = runner.invoke(app, ["split", str(f), "-f", "train=0.7,test=0.3", "-s", "identity", "--seed", "3", "--endpoint-column", "y", "-o", str(out), "--json", str(js)])
    assert r.exit_code == 0, r.output
    assert "indivisible groups" in r.output and "verification train→test" in r.output
    rows = list(csv.DictReader(open(out)))
    assert len(rows) == len(SMILES) and {r_["split"] for r_ in rows} == {"train", "test"}
    assert json.loads(js.read_text())["leakage"]["identity_level"] == "parent"
    assert runner.invoke(app, ["split", str(f), "-s", "nope"]).exit_code == 1
    assert runner.invoke(app, ["split", str(f), "-f", "train=0.5,test=0.4"]).exit_code == 1
    assert runner.invoke(app, ["split", str(f), "-f", "bad"]).exit_code == 1
    assert runner.invoke(app, ["split", str(f), "-s", "temporal"]).exit_code == 1
    assert runner.invoke(app, ["split", str(f), "-s", "temporal", "--date-column", "date"]).exit_code == 0

    gen = tmp_path / "gen.smi"; gen.write_text("CCO\nCCO\nbad((\n\nc1ccccc1\n")
    ref = tmp_path / "train.smi"; ref.write_text("CCO\n")
    gout, gjs = tmp_path / "g.csv", tmp_path / "g.json"
    r = runner.invoke(app, ["generated", str(gen), "-r", str(ref), "--constraints", "mw=0:200", "--repair", "-o", str(gout), "--json", str(gjs)])
    assert r.exit_code == 0, r.output
    assert "Validity" in r.output and "unique valid" in r.output and "attempts" in r.output
    data = json.loads(gjs.read_text())
    assert data["n_generated"] == 5 and data["n_valid"] == 3 and data["reference_sets"] == ["train"]
    assert len(list(csv.DictReader(open(gout)))) == 5
    assert runner.invoke(app, ["generated", str(gen), "--constraints", "oops"]).exit_code == 1
