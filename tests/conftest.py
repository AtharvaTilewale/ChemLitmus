"""
Pytest configuration and shared fixtures.

Provides fixtures for testing, including temporary directories, mock clients, etc.
"""

import tempfile
from pathlib import Path

import pytest

from chemlitmus.config import Settings


@pytest.fixture
def temp_cache_dir() -> Path:
    """Provide a temporary cache directory for testing."""
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def test_settings(temp_cache_dir: Path) -> Settings:
    """Provide test settings with temporary directories."""
    return Settings(
        cache_dir=temp_cache_dir,
        data_dir=temp_cache_dir / "data",
        log_dir=temp_cache_dir / "logs",
    )


@pytest.fixture
def mock_smiles_list() -> list[str]:
    """Provide a list of valid SMILES for testing."""
    return [
        "c1ccccc1",  # Benzene
        "CCO",  # Ethanol
        "CC(=O)O",  # Acetic acid
        "c1ccc(O)cc1",  # Phenol
    ]
