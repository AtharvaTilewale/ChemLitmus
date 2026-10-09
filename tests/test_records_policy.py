"""Record-preserving ingestion, policy hashing and manifests."""

import json

import pytest
from rdkit import Chem

from chemlitmus import ChemicalPolicy, SchemaError, read_records, resolve_policy
from chemlitmus.core.manifest import new_manifest, sha256_text
from chemlitmus.core.policy import RepairPolicy
from chemlitmus.utils.parsers import parse_compounds_file


def test_table_accounting_and_metadata(tmp_path):
    f = tmp_path / "d.csv"
    f.write_text("id,smiles,activity\nA1,CCO,1.2\nA2,,3\nA3,CC O,4\nA4,C((,5\nA1,c1ccccc1,6\n")
    rs = read_records(f)
    assert rs.roles == {"structure": "smiles", "id": "id", "endpoint": "activity"}
    assert rs.n_total == 5 and rs.n_ok == 2 and rs.n_empty == 1 and rs.n_invalid == 2 and rs.reconcile()
    r = {x.position: x for x in rs.records}
    assert r[0].source_id == "A1" and r[0].fields == {"id": "A1", "activity": "1.2"} and r[0].parsed_smiles == "CCO"
    assert r[1].status == "empty" and r[1].issues[0].code == "PARSE_EMPTY"
    assert r[2].status == "invalid" and r[2].issues[0].code == "PARSE_WHITESPACE" and r[2].issues[0].evidence["first_token"] == "CC"
    assert r[3].status == "invalid" and r[3].issues[0].code == "PARSE_INVALID"
    assert r[4].source_id == "A1" and r[4].record_id != r[0].record_id        # duplicate source ids stay distinct
    assert [x.record_id for x in rs.records] == [f"r{i:06d}" for i in range(5)]


def test_ambiguous_columns_are_a_schema_error(tmp_path):
    f = tmp_path / "d.csv"; f.write_text("foo,bar\nCCO,x\n")
    with pytest.raises(SchemaError):
        read_records(f)
    rs = read_records(f, ambiguous="first")
    assert rs.roles["structure"] == "foo" and rs.notes
    rs = read_records(f, structure_column="foo", id_column="bar")
    assert rs.records[0].source_id == "x"
    with pytest.raises(SchemaError):
        read_records(f, structure_column="nope")
    # header-less single column is unambiguous; header-less multi-column is not
    (tmp_path / "s.csv").write_text("CCO\nCCN\n")
    assert read_records(tmp_path / "s.csv").n_ok == 2
    (tmp_path / "m.csv").write_text("CCO,ethanol\nCCN,ethylamine\n")
    with pytest.raises(SchemaError):
        read_records(tmp_path / "m.csv")


def test_explicit_roles(tmp_path):
    f = tmp_path / "d.tsv"; f.write_text("cpd\tval\tu\nCCO\t5\tnM\n")
    rs = read_records(f, structure_column="cpd", roles={"endpoint": "val", "units": "u"})
    assert rs.roles == {"structure": "cpd", "endpoint": "val", "units": "u"} and rs.format == "tsv"
    with pytest.raises(SchemaError):
        read_records(f, structure_column="cpd", roles={"endpoint": "missing"})


def test_smi_keeps_names_and_comments(tmp_path):
    f = tmp_path / "x.smi"; f.write_text("# header comment\nCCO ethanol\nCCN\n\nC(( broken\n")
    rs = read_records(f)
    s = {x.position: x for x in rs.records}
    assert s[0].status == "unsupported" and s[1].source_id == "ethanol" and s[1].structure == "CCO"
    assert s[2].source_id is None and s[3].status == "empty" and s[4].status == "invalid" and s[4].source_id == "broken"
    assert rs.n_total == 5 and rs.reconcile()


def test_sdf_keeps_properties_and_failed_records(tmp_path):
    good = Chem.MolToMolBlock(Chem.MolFromSmiles("CCO"))
    good = "ethanol" + good[good.index("\n"):]                                     # title line
    bad = good.replace("  3  2  0  0  0  0  0  0  0  0999 V2000", "  3  9  0  0  0  0  0  0  0  0999 V2000")  # corrupt bond count
    text = good + ">  <ID>\nE1\n\n>  <Activity>\n4.2\n\n$$$$\n" + bad.replace("ethanol", "broken") + ">  <ID>\nE2\n\n$$$$\n"
    f = tmp_path / "m.sdf"; f.write_text(text)
    rs = read_records(f)
    assert rs.n_total == 2 and rs.n_ok == 1 and rs.n_invalid == 1 and rs.reconcile()
    ok, broken = rs.records
    assert ok.source_id == "ethanol" and ok.fields == {"ID": "E1", "Activity": "4.2"} and ok.parsed_smiles == "CCO" and ok.structure_format == "molblock"
    assert broken.status == "invalid" and broken.issues[0].code == "SDF_PARSE_FAILED" and broken.source_id == "broken" and broken.fields == {"ID": "E2"}
    assert rs.columns == ["title", "ID", "Activity"]
    assert parse_compounds_file(f) == ["CCO"]                                     # legacy wrapper: parsed molecules only


def test_xlsx(tmp_path):
    import pandas as pd
    f = tmp_path / "d.xlsx"
    pd.DataFrame({"Name": ["a", "b"], "SMILES": ["CCO", "bad(("]}).to_excel(f, index=False)
    rs = read_records(f)
    assert rs.roles == {"structure": "SMILES", "id": "Name"} and rs.n_ok == 1 and rs.n_invalid == 1


def test_policy_presets_hash_and_roundtrip(tmp_path):
    p = ChemicalPolicy.preset("parent"); c = ChemicalPolicy.preset("conservative")
    assert p.hash != c.hash and p.hash == ChemicalPolicy.preset("parent").hash and len(p.hash) == 64
    assert c.fragment == "keep_all" and c.identity_level == "exact" and c.repair.mode == "none"
    f = tmp_path / "pol.json"; p.save(f)
    assert ChemicalPolicy.load(f).hash == p.hash and resolve_policy(f).hash == p.hash and resolve_policy("conservative").hash == c.hash
    assert resolve_policy(None).name == "parent" and resolve_policy({"identity_level": "skeleton"}).identity_level == "skeleton"
    with pytest.raises(ValueError):
        ChemicalPolicy(identity_level="formula")
    with pytest.raises(ValueError):
        ChemicalPolicy(repair=RepairPolicy(mode="maybe"))
    with pytest.raises(ValueError):
        ChemicalPolicy.preset("nope")


def test_manifest(tmp_path):
    f = tmp_path / "in.csv"; f.write_text("smiles\nCCO\n")
    m = new_manifest("chemlitmus audit", structure_column="smiles")
    m.add_input(f); m.policy_hash = ChemicalPolicy.preset("parent").hash
    out = tmp_path / "manifest.json"; m.save(out)
    data = json.loads(out.read_text())
    assert data["inputs"][0]["path"] == "in.csv" and len(data["inputs"][0]["sha256"]) == 64 and data["rdkit_version"]
    assert "/" not in data["inputs"][0]["path"] and data["settings"] == {"structure_column": "smiles"}
    assert sha256_text("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
