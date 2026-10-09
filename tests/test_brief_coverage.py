"""Curated regression corpus and the verification checks the brief asks for.

Every expectation below states *why* it holds. Structures are either textbook cases (ethanol,
benzene, glycine) or public reference compounds named in the comment; no proprietary data is used.
"""

import csv
import json

import pytest
from rdkit import Chem

from chemlitmus import ChemicalPolicy, audit_dataset, evaluate_generated, label_conflicts, make_splits, read_records
from chemlitmus.core.cleanup import propose_cleanup
from chemlitmus.core.identity import compute_identity
from chemlitmus.core.smartsaudit import audit_smarts
from chemlitmus.core.standardize import standardize_smiles

# --------------------------------------------------------------------------- curated corpus
# (label, structure, expected status, reason)
CORPUS = [
    ("valid", "CCO", "ok", "ethanol: textbook valid SMILES"),
    ("blank", "", "empty", "blank cell must be preserved as a record, not dropped"),
    ("whitespace", "CC O", "invalid", "RDKit would parse ethane and treat 'O' as a title"),
    ("truncated_ring", "C1CCCCC", "invalid", "unclosed ring bond 1"),
    ("truncated_branch", "CC(C", "invalid", "unbalanced parenthesis"),
    ("bad_valence", "CN(C)(C)C", "invalid", "five-coordinate neutral nitrogen is rejected by sanitisation"),
    ("salt", "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]", "ok", "aspirin sodium: two components, net neutral"),
    ("mixture", "O=C(O)c1ccccc1O.CC(=O)Oc1ccccc1C(=O)O", "ok", "salicylic acid + aspirin: two substantial organic components"),
    ("zwitterion", "[NH3+]CC(=O)[O-]", "ok", "glycine zwitterion: charges balance, net neutral"),
    ("permanent_charge", "C[N+](C)(C)C", "ok", "tetramethylammonium: permanent cation, no neutral form"),
    ("isotope", "[13CH3]CO", "ok", "carbon-13 labelled ethanol"),
    ("tautomer_keto", "O=C1CCCCC1", "ok", "cyclohexanone, keto form"),
    ("tautomer_enol", "OC1=CCCCC1", "ok", "cyclohexenol, enol form of the same compound"),
    ("stereo_assigned", "C[C@H](N)C(=O)O", "ok", "L-alanine, stereocentre assigned"),
    ("stereo_unassigned", "CC(N)C(=O)O", "ok", "alanine with the centre undefined"),
    ("stereo_partial", "C[C@H](O)C(O)CO", "ok", "one centre assigned, one not"),
    ("acyclic", "CCCCCC", "ok", "hexane: no Murcko scaffold"),
    ("organometallic", "[Fe+2].[C-]#[O+]", "ok", "iron carbonyl fragment: parses, but outside descriptor/alert scope"),
]


@pytest.fixture
def corpus_csv(tmp_path):
    f = tmp_path / "corpus.csv"
    with open(f, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "smiles", "note"])
        for label, smi, _status, reason in CORPUS:
            w.writerow([label, smi, reason])          # quoted fields with commas exercise the CSV reader
    return f


def test_corpus_status_and_accounting(corpus_csv):
    rs = read_records(corpus_csv)
    got = {r.source_id: r.status for r in rs.records}
    for label, _smi, expected, reason in CORPUS:
        assert got[label] == expected, f"{label}: expected {expected} ({reason}), got {got[label]}"
    assert rs.reconcile() and rs.n_total == len(CORPUS)
    assert rs.records[0].fields["note"].startswith("ethanol")      # quoted field with a comma survived


def test_corpus_audit_flags(corpus_csv):
    a = audit_dataset(corpus_csv)
    codes = {r.source_id: {i.code for i in r.issues} for r in a.records}
    assert "FRAG_MULTIPLE_ORGANIC" in codes["mixture"] and "FRAG_MULTIPLE_ORGANIC" not in codes["salt"]
    assert "FRAG_MULTIPLE" in codes["salt"] and "FRAG_MULTIPLE" not in codes["zwitterion"]   # a zwitterion is one component
    assert "CHARGE_NET" in codes["permanent_charge"] and "CHARGE_NET" not in codes["zwitterion"]
    assert "ISOTOPE" in codes["isotope"]
    assert "STEREO_UNASSIGNED" in codes["stereo_unassigned"] and "STEREO_PARTIAL" in codes["stereo_partial"]
    assert "STEREO_UNASSIGNED" not in codes["stereo_assigned"] and "STEREO_PARTIAL" not in codes["stereo_assigned"]
    assert "ELEMENT_UNUSUAL" in codes["organometallic"]
    assert "PARSE_WHITESPACE" in codes["whitespace"] and "PARSE_EMPTY" in codes["blank"]
    assert {"PARSE_INVALID"} <= codes["truncated_ring"] & codes["truncated_branch"] & codes["bad_valence"]
    assert a.summary.n_total == len(CORPUS) and a.summary.n_annotated == sum(1 for *_x, s, _r in [(c[1], c[2], c[3]) for c in CORPUS] if s == "ok")
    by_id = {r.source_id: r for r in a.records}
    assert by_id["acyclic"].scaffold == "acyclic"


# --------------------------------------------------------------------------- invariants the brief asks for

@pytest.mark.parametrize("smi", [c[1] for c in CORPUS if c[2] == "ok"])
def test_standardization_is_idempotent(smi):
    """Standardising twice under the same policy must equal standardising once."""
    once = standardize_smiles(smi)
    assert once.success
    twice = standardize_smiles(once.output_smiles)
    assert twice.success and twice.output_smiles == once.output_smiles


@pytest.mark.parametrize("a,b", [
    ("CCO", "OCC"), ("c1ccccc1", "C1=CC=CC=C1"), ("CC(=O)Oc1ccccc1C(=O)O", "O=C(O)c1ccccc1OC(C)=O"),
    ("C[C@H](N)C(=O)O", "N[C@@H](C)C(=O)O"), ("[NH3+]CC(=O)[O-]", "[O-]C(=O)C[NH3+]"),
])
def test_identity_is_invariant_under_equivalent_smiles(a, b):
    """Two spellings of one molecule must give identical keys at every level."""
    ka, kb = compute_identity(a), compute_identity(b)
    assert ka.is_valid and kb.is_valid
    for level in ("exact", "parent", "tautomer", "nostereo", "skeleton", "formula"):
        assert ka.key(level) == kb.key(level), level


def test_splits_preserve_groups_and_seed(corpus_csv):
    rs = read_records(corpus_csv)
    recs = [{"record_id": r.record_id, "smiles": r.parsed_smiles or "", "source_id": r.source_id} for r in rs.records]
    r1 = make_splits(recs, {"train": 0.6, "test": 0.4}, strategy="identity", seed=5)
    r2 = make_splits(recs, {"train": 0.6, "test": 0.4}, strategy="identity", seed=5)
    assert [a.split for a in r1.assignments] == [a.split for a in r2.assignments]
    by_group = {}
    for a in r1.assignments:
        by_group.setdefault(a.group_key, set()).add(a.split)
    assert all(len(v) == 1 for v in by_group.values()), "a group was divided"


def test_cli_and_api_agree(corpus_csv, tmp_path):
    """The same input must give the same chemical decisions through both interfaces."""
    from typer.testing import CliRunner

    from chemlitmus.cli import app

    out = tmp_path / "o"
    r = CliRunner().invoke(app, ["audit", str(corpus_csv), "-o", str(out)])
    assert r.exit_code == 0, r.output
    cli_summary = json.loads((out / "summary.json").read_text())
    api = audit_dataset(corpus_csv)
    assert cli_summary["n_total"] == api.summary.n_total and cli_summary["n_annotated"] == api.summary.n_annotated
    assert cli_summary["issues_by_code"] == api.summary.issues_by_code
    assert cli_summary["policy_hash"] == api.policy_hash and cli_summary["n_groups"] == api.summary.n_groups


# --------------------------------------------------------------------------- newly added brief items

def test_technical_vs_biological_replicates():
    """§8: repeated entries from one source are technical; disagreement across sources is not."""
    recs = [
        {"record_id": "1", "key": "K", "val": "10", "u": "nM", "doc": "D1"},
        {"record_id": "2", "key": "K", "val": "12", "u": "nM", "doc": "D1"},    # same experiment
        {"record_id": "3", "key": "K", "val": "5000", "u": "nM", "doc": "D2"},  # different experiment
    ]
    rep = label_conflicts(recs, endpoint_field="val", units_field="u", source_field="doc")
    g = rep.groups[0]
    assert rep.source_field == "doc" and g.conflict
    rp = g.replicates
    assert rp.basis == "source field" and rp.n_sources == 2 and rp.n_technical_groups == 1 and rp.n_technical_records == 2
    assert rp.technical_spread < 0.1 and rp.between_source_spread > 2    # ~0.08 log within D1, ~2.6 log between
    no_src = label_conflicts(recs, endpoint_field="val", units_field="u")
    assert no_src.groups[0].replicates.basis == "not determinable"


def test_leakage_states_its_scope():
    """§7: measured contamination against supplied collections is not evidence about pretraining."""
    from chemlitmus import leakage_report

    rep = leakage_report({"train": [("1", "CCO", None)], "test": [("2", "CCO", None)]})
    assert "pretraining" in rep.scope and "supplied" in rep.scope and "cannot see" in rep.scope


def test_cleanup_proposals_preserve_provenance():
    """§12: redundancy across published sets is provenance, never an automatic removal."""
    pats = [("[OX2H]", "alcohol", "SetA"), ("[OX2H]", "alcohol-copy", "SetA"), ("[CX4][OX2H]", "alkylol", "SetB")]
    lib = [Chem.MolFromSmiles(s) for s in ["CCO", "c1ccccc1O", "CC(=O)O"]]
    res = audit_smarts(pats, library=lib, checks=["compile", "breadth", "dead", "redundancy", "proof"])
    rep = propose_cleanup(res)
    by_name = {p.name: p for p in rep.proposals}
    # one member of the equivalent pair is kept; only the duplicate is proposed
    assert "alcohol" not in by_name and rep.n_retained_representatives >= 1
    assert by_name["alcohol-copy"].action == "recommended" and "exact text duplicate" in by_name["alcohol-copy"].evidence
    assert by_name["alkylol"].action == "keep provenance" and by_name["alkylol"].crosses_rule_sets
    assert "SetA" in by_name["alkylol"].covered_by_rule_sets and "provenance" in rep.note
    assert propose_cleanup(res, allow_cross_set=True).by_action.get("keep provenance") is None
    assert rep.n_verdict_changes_if_all_applied == 0


def test_provider_records_carry_retrieval_provenance():
    """§13: a merged record must say which source each field came from, and when it was fetched."""
    from chemlitmus.providers import CompoundRecord, _merge

    a = CompoundRecord(source="chembl", source_id="CHEMBL25", query="aspirin", name="ASPIRIN", formula="C9H8O4", retrieved_at="2026-10-09T10:00:00+00:00")
    b = CompoundRecord(source="chebi", source_id="CHEBI:15365", query="aspirin", smiles="CC(=O)Oc1ccccc1C(=O)O", xlogp=1.2,
                       retrieved_at="2026-10-09T11:00:00+00:00", from_cache=True, cached_at="2026-10-01T00:00:00+00:00")
    m = _merge([a, b], "aspirin")
    assert m.field_provenance["name"] == "chembl" and m.field_provenance["smiles"] == "chebi" and m.field_provenance["xlogp"] == "chebi"
    assert m.retrieved_at == "2026-10-09T11:00:00+00:00" and m.from_cache is False      # not every part came from cache
    assert b.cached_at != b.retrieved_at                                                  # stale cache is distinguishable


def test_parquet_roundtrip(tmp_path):
    """§5: Parquet is an optional adapter with the same record semantics."""
    pytest.importorskip("pyarrow")
    import pandas as pd

    f = tmp_path / "d.parquet"
    pd.DataFrame({"id": ["A", "B", "C"], "smiles": ["CCO", "bad((", ""]}).to_parquet(f)
    rs = read_records(f)
    assert rs.format == "parquet" and rs.roles == {"structure": "smiles", "id": "id"}
    assert (rs.n_ok, rs.n_invalid, rs.n_empty) == (1, 1, 1) and rs.reconcile()
    assert rs.records[0].parsed_smiles == "CCO" and rs.records[0].source_id == "A"


def test_generation_reports_reference_scope():
    g = evaluate_generated(["CCO"], reference={"train": ["CCN"]})
    assert g.reference_sets == ["train"] and "reference sets listed" in g.note


def test_conservative_policy_keeps_everything(corpus_csv):
    a = audit_dataset(corpus_csv, policy=ChemicalPolicy.preset("conservative"))
    by_id = {r.source_id: r for r in a.records}
    assert by_id["salt"].standardized_smiles == by_id["salt"].parsed_smiles      # no fragment choice, no neutralisation
    assert by_id["zwitterion"].standardized_smiles == by_id["zwitterion"].parsed_smiles
