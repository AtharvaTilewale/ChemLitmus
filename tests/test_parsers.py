"""Regression tests for input file parsing (header detection must not eat data rows)."""

from pathlib import Path

import pytest

from chemlitmus.utils.parsers import parse_compounds_file


@pytest.mark.parametrize("name,content,expected", [
    ("names.txt", "aspirin\nCHEMBL25\n", ["aspirin", "CHEMBL25"]),
    ("ids.txt", "2244\n5090\n", ["2244", "5090"]),
    ("smiles.txt", "CCO\nCCN\n", ["CCO", "CCN"]),
    ("header.txt", "smiles\nCCO\n", ["CCO"]),
    ("bom.txt", "\ufeffCCO\nCCN\n", ["CCO", "CCN"]),
    ("multi.csv", "id,smiles,name\n1,CCO,ethanol\n2,c1ccccc1,benzene\n", ["CCO", "c1ccccc1"]),
    ("headerless_multi.csv", "CCO,ethanol\nCCN,ethylamine\n", ["CCO", "CCN"]),
    ("names_header.csv", "Name\naspirin\ncaffeine\n", ["aspirin", "caffeine"]),
    ("cid_header.csv", "CID\n2244\n5090\n", ["2244", "5090"]),
    ("tsv.tsv", "name\tcanonical_smiles\nx\tCCO\n", ["CCO"]),
    ("blank_lines.txt", "CCO\n\nCCN\n", ["CCO", "CCN"]),
    ("smi.smi", "CCO ethanol\nCCN ethylamine\n", ["CCO", "CCN"]),
])
def test_parse_compounds_file(tmp_path: Path, name, content, expected):
    f = tmp_path / name
    f.write_text(content, encoding="utf-8")
    assert parse_compounds_file(f) == expected


def test_missing_and_unsupported(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        parse_compounds_file(tmp_path / "nope.csv")
    f = tmp_path / "x.xyz"; f.write_text("CCO")
    with pytest.raises(ValueError):
        parse_compounds_file(f)
