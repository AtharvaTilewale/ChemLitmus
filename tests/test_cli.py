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
