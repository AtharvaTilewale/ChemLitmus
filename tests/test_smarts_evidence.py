"""Evidence labels, match semantics, catalogue metadata, reference-panel summary and holdout."""

from rdkit import Chem
from typer.testing import CliRunner

from chemlitmus.cli import app
from chemlitmus.core.smartsaudit import EVIDENCE_LABELS, CatalogueMetadata, audit_smarts

runner = CliRunner()
LIB = [Chem.MolFromSmiles(s) for s in ["CCO", "c1ccccc1", "CC(=O)O", "c1ccncc1", "CCN", "C1CCCCC1"]]
HOLD = [Chem.MolFromSmiles(s) for s in ["ClCCl", "c1ccc(Cl)cc1"]]
PATS = [("[OX2H]", "alcohol", "A"), ("[OX2H]", "dup", "A"), ("[Cl]", "chloro", "B"), ("[CX4][OX2H]", "alkylol", "B")]


def test_evidence_labels_distinguish_observed_from_proven():
    res = audit_smarts(PATS, library=LIB, library_source="demo", checks=["compile", "breadth", "dead", "redundancy", "sensitivity", "proof"])
    ev = {p.name: p.evidence for p in res.patterns}
    assert ev["dup"][0] == "exact text duplicate" and "identical observed hit set" in ev["dup"] and "static proven equivalence" in ev["dup"]
    assert ev["alkylol"] == ["observed hit-set containment", "static proven containment"]
    assert ev["chloro"] == ["undecided", "not observed in this reference"]        # never 'impossible'
    for p in res.patterns:
        assert set(p.evidence) <= set(EVIDENCE_LABELS)


def test_examples_and_matched_atoms():
    res = audit_smarts(PATS, library=LIB, library_source="demo", n_examples=2)
    a = res.patterns[0]
    assert a.example_matches == ["CCO", "CC(=O)O"] and a.example_match_atoms == [[2], [3]]
    assert res.patterns[2].example_matches == []
    rows = res.to_rows()
    assert rows[0]["example_matches"].startswith("CCO") and "observed" in rows[3]["evidence"] or rows[3]["evidence"] == ""


def test_match_semantics_catalogue_and_panel_are_recorded():
    res = audit_smarts(PATS, library=LIB, library_source="demo", catalogue=CatalogueMetadata(name="demo", version="1", licence="CC0", source="unit test"))
    ms = res.match_semantics
    assert ms.preparation == "implicit-h" and ms.use_chirality is False and ms.hydrogens == "implicit" and ms.rdkit_version
    c = res.catalogue
    assert c.name == "demo" and c.licence == "CC0" and c.n_patterns == 4 and c.rule_sets == {"A": 2, "B": 2} and len(c.content_hash) == 64
    p = res.reference_panel
    assert p.n_molecules == 6 and p.median_heavy_atoms == 6.0 and p.fraction_with_ring == 0.5 and p.elements["C"] == 6 and len(p.content_hash) == 64
    # the same patterns in a different order hash the same (content, not order)
    assert audit_smarts(list(reversed(PATS)), library=LIB).catalogue.content_hash == c.content_hash


def test_holdout_revives_patterns_not_observed_on_the_main_panel():
    res = audit_smarts(PATS, library=LIB, library_source="demo", holdout_libraries={"halogens": HOLD})
    h = res.holdout[0]
    assert h.panel == "halogens" and h.n_molecules == 2 and h.n_patterns_firing == 1
    assert h.patterns_revived == [2] and h.hits == {2: 2}
    assert res.patterns[2].never_fires is True                                   # still 'not observed here' — the holdout is the evidence against extrapolating
    assert "not observed in this reference" in res.patterns[2].evidence


def test_preparation_flip_examples():
    pats = [("[#6][H]", "needs explicit H", None)]
    res = audit_smarts(pats, library=LIB, library_source="demo", checks=["compile", "breadth", "dead", "sensitivity"])
    p = res.patterns[0]
    assert p.n_hits == 0 and p.dead_verdict and "explicit-h" in p.dead_verdict
    assert p.preparation_flip_examples.get("explicit-h")


def test_cli_passes_catalogue_and_holdout(tmp_path):
    pf = tmp_path / "p.csv"; pf.write_text("name,smarts,rule_set\nalcohol,[OX2H],A\nchloro,[Cl],B\n")
    lf = tmp_path / "lib.smi"; lf.write_text("CCO\nc1ccccc1\n")
    hf = tmp_path / "hold.smi"; hf.write_text("ClCCl\n")
    r = runner.invoke(app, ["smartsaudit", str(pf), "--library", str(lf), "--holdout", str(hf), "--catalogue-name", "demo", "--catalogue-licence", "CC0"])
    assert r.exit_code == 0, r.output
    assert "Holdout hold" in r.output and "panel-specific" in r.output
    assert "Catalogue: demo" in r.output and "not one named set" in r.output
    assert "Match semantics" in r.output and "Reference panel" in r.output
    assert runner.invoke(app, ["smartsaudit", str(pf), "--holdout", str(tmp_path / "nope.smi")]).exit_code == 1
