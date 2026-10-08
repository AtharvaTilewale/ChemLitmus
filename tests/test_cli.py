"""Tests for the Typer CLI application."""

import re
from pathlib import Path
import pytest
from typer.testing import CliRunner
from chemlitmus.cli.main import app

# Initialize without the invalid mix_stderr argument
runner = CliRunner()


def _clean(text: str) -> str:
    """Strip ANSI terminal escape codes for reliable string assertions."""
    return re.sub(r"\x1b\[[0-9;]*[a-zA-Z]", "", text)


def test_app_version():
    """Test the version flag."""
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "ChemLitmus" in clean_out
    assert "License:" in clean_out
    assert "MIT" in clean_out
    assert "Repository:" in clean_out
    assert "AtharvaTilewale/ChemLitmus" in clean_out


def test_app_status():
    """Test the status command."""
    result = runner.invoke(app, ["status"])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "ChemLitmus Status" in clean_out
    assert "Configuration" in clean_out
    assert "Cache Dir" in clean_out


def test_app_lookup_missing_argument():
    """Test lookup command fails gracefully when missing query."""
    result = runner.invoke(app, ["lookup"])
    assert result.exit_code != 0
    clean_out = _clean(result.output)
    assert "Missing argument" in clean_out


@pytest.mark.integration
def test_app_lookup_compound_name():
    """Test lookup command for compound name returns SMILES (requires network)."""
    result = runner.invoke(app, ["lookup", "aspirin"])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "CC(=O)OC1=CC=CC=C1C(=O)O" in clean_out
    assert "2244" in clean_out


def test_download_help_includes_gen():
    """Test that download --help exposes the --gen option."""
    result = runner.invoke(app, ["download", "--help"])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "--gen" in clean_out or "-g" in clean_out


def test_download_gen_all_single_smiles(tmp_path: Path):
    """Test download with --gen all for a single SMILES in 2D and 3D."""
    out_dir = tmp_path / "structs"

    # 3D SDF
    result = runner.invoke(app, ["download", "CCO", "--gen", "all", "--3d", "--format", "sdf", "-o", str(out_dir)])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "Successfully generated" in clean_out

    # 2D MOL
    result = runner.invoke(app, ["download", "c1ccccc1", "--gen", "all", "--2d", "--format", "mol", "-o", str(out_dir)])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "Successfully generated" in clean_out

    # 3D PDB
    result = runner.invoke(app, ["download", "CC(=O)O", "--gen", "all", "--3d", "--format", "pdb", "-o", str(out_dir)])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "Successfully generated" in clean_out


def test_download_gen_all_batch_file(tmp_path: Path):
    """Test download with --gen all on a batch file of SMILES."""
    smi_file = tmp_path / "compounds.smi"
    smi_file.write_text("c1ccccc1\nCCO\nCC(=O)O\n", encoding="utf-8")
    out_dir = tmp_path / "batch_out"

    result = runner.invoke(app, ["download", "--file", str(smi_file), "--gen", "all", "--3d", "--format", "sdf", "-o", str(out_dir)])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "Batch processing complete" in clean_out
    assert "Generated:" in clean_out


@pytest.mark.integration
def test_download_gen_missing_batch_file(tmp_path: Path):
    """Test download with --gen missing generates structures not found in PubChem (requires network)."""
    smi_file = tmp_path / "compounds_missing.smi"
    smi_file.write_text("c1ccccc1\nCC(C)(C)CC(=O)N1CCCCC1C(=O)O\n", encoding="utf-8")
    out_dir = tmp_path / "batch_missing_out"

    result = runner.invoke(app, ["download", "--file", str(smi_file), "--gen", "missing", "--3d", "--format", "sdf", "-o", str(out_dir)])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "Batch processing complete" in clean_out


def test_download_gen_invalid_mode():
    """Test invalid --gen argument."""
    result = runner.invoke(app, ["download", "CCO", "--gen", "invalid_mode"])
    assert result.exit_code != 0
    clean_out = _clean(result.output)
    assert "Invalid value for --gen" in clean_out


def test_app_update_help():
    """Test update --help."""
    result = runner.invoke(app, ["update", "--help"])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "Check for a newer version" in clean_out
    assert "--check" in clean_out or "-c" in clean_out
    assert "--yes" in clean_out or "-y" in clean_out


def test_detect_install_source():
    """Test install source detector helper."""
    from chemlitmus.cli.main import _detect_install_source
    source_type, detail = _detect_install_source()
    assert source_type in ["git_repo", "git_pip", "pip"]
    assert detail is not None


def test_fingerprint_single_cli():
    """Test fingerprint command for a single SMILES."""
    result = runner.invoke(app, ["fingerprint", "CC(=O)OC1=CC=CC=C1C(=O)O", "--type", "ecfp4"])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "ECFP4" in clean_out
    assert "2048" in clean_out


def test_filter_single_cli():
    """Test filter command for a single SMILES."""
    result = runner.invoke(app, ["filter", "CC(=O)OC1=CC=CC=C1C(=O)O", "--rules", "lipinski,veber"])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "Lipinski Ro5" in clean_out
    assert "Veber" in clean_out


def test_similar_cli(tmp_path: Path):
    """Test similar command with a library file."""
    lib_file = tmp_path / "lib.smi"
    lib_file.write_text("CCO\nCCCO\nCC(=O)OC1=CC=CC=C1C(=O)O\n", encoding="utf-8")
    result = runner.invoke(app, ["similar", "CC(=O)OC1=CC=CC=C1C(=O)O", "--file", str(lib_file), "--threshold", "0.1"])
    assert result.exit_code == 0
    clean_out = _clean(result.output)
    assert "Similarity Search Results" in clean_out


# --------------------------------------------------------------------------- regression tests (review pass)

def test_validate_command_single_and_quiet():
    r = runner.invoke(app, ["validate", "CCO"])
    assert r.exit_code == 0 and "C2H6O" in r.output
    r = runner.invoke(app, ["validate", "CCO", "--quiet"])
    assert r.exit_code == 0 and r.output.strip() == "CCO"
    r = runner.invoke(app, ["validate", "C1CC(C"])
    assert r.exit_code == 2 and "diagnose" in r.output
    assert runner.invoke(app, ["validate"]).exit_code == 1


def test_validate_command_batch(tmp_path):
    f = tmp_path / "lib.smi"; f.write_text("CCO\nC1CC(C\nc1ccccc1\n")
    out = tmp_path / "val.csv"
    r = runner.invoke(app, ["validate", "--file", str(f), "--output", str(out)])
    assert r.exit_code == 0 and "valid=2" in r.output and "invalid=1" in r.output
    assert out.read_text().count("\n") == 4


def test_tautomers_invalid_smiles_exits_nonzero():
    assert runner.invoke(app, ["tautomers", "C1CC(C"]).exit_code == 1


def test_standardize_bad_step_single_clean_error():
    r = runner.invoke(app, ["standardize", "CCO", "--steps", "bogus"])
    assert r.exit_code == 1 and "Invalid steps" in r.output
    assert "Error: \n" not in r.output


def test_filter_bad_prep_is_usage_error(tmp_path):
    f = tmp_path / "lib.smi"; f.write_text("CCO\n")
    r = runner.invoke(app, ["filter", "--file", str(f), "--rules", "lipinski", "--prep", "bogus"])
    assert r.exit_code == 1 and "Unknown preparation" in r.output


def test_lookup_bad_type_rejected_before_network():
    r = runner.invoke(app, ["lookup", "CCO", "--type", "bogus"])
    assert r.exit_code == 1 and "Unknown --type" in r.output


def test_conformers_num_alias_and_default_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    r = runner.invoke(app, ["conformers", "CCO", "--num", "2"])
    assert r.exit_code == 0 and (tmp_path / "conformers.sdf").exists()


def test_rgroup_positional_core(tmp_path):
    f = tmp_path / "lib.smi"; f.write_text("Cc1ccccc1\nCCc1ccccc1\n")
    r = runner.invoke(app, ["rgroup", "c1ccccc1[*:1]", "--file", str(f)])
    assert r.exit_code == 0 and "R-Group" in r.output
    assert runner.invoke(app, ["rgroup", "--file", str(f)]).exit_code == 1


def test_scaffold_acyclic_is_annotated():
    r = runner.invoke(app, ["scaffold", "CCO"])
    assert r.exit_code == 0 and "acyclic" in r.output


def test_download_gen_all_is_offline(tmp_path, monkeypatch):
    import chemlitmus.cli.main as m
    def boom(*a, **k): raise AssertionError("network lookup attempted in --gen all")
    monkeypatch.setattr(m, "lookup", boom)
    r = runner.invoke(app, ["download", "CCO", "--gen", "all", "--2d", "--format", "mol", "--output-dir", str(tmp_path)])
    assert r.exit_code == 0, r.output
    assert any(p.suffix == ".mol" for p in tmp_path.iterdir())


def test_diagnose_batch_does_not_leak_rdkit_log(tmp_path):
    f = tmp_path / "lib.smi"; f.write_text("CCO\nC1CC(C\nCN(C)(C)C\n")
    r = runner.invoke(app, ["diagnose", "--file", str(f)])
    assert r.exit_code == 0
    assert "SMILES Parse Error" not in r.output
