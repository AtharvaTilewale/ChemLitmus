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


@pytest.fixture(autouse=True)
def _isolated_cache(tmp_path, monkeypatch):
    """Every test gets its own SQLite cache so the suite never touches the user's cache."""
    from chemlitmus.config import settings
    import chemlitmus.providers.base as pb
    cache = tmp_path / "cache"
    cache.mkdir()
    monkeypatch.setattr(settings, "cache_dir", cache)
    monkeypatch.setattr(pb, "_DB", None)
    yield
    monkeypatch.setattr(pb, "_DB", None)
