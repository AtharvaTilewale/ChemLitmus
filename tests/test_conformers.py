"""Tests for chemlitmus.core.conformers"""
import os
from chemlitmus.core.conformers import generate_conformers

def test_generate_conformers_valid(tmp_path):
    sdf_out = str(tmp_path / "out.sdf")
    result = generate_conformers("CCO", num_conformers=3, output_sdf=sdf_out)
    assert result.error is None
    assert result.num_generated > 0
    assert result.num_generated <= 3
    assert os.path.exists(sdf_out)

def test_generate_conformers_invalid():
    result = generate_conformers("invalid")
    assert result.error is not None
    assert result.num_generated == 0
