"""Tests for SQLite caching mechanism."""

from pathlib import Path
from chemlitmus.core.database import DatabaseManager
from chemlitmus.core.pubchem import PubChemCompound


def test_database_init(temp_cache_dir: Path) -> None:
    """Test SQLite database initialization."""
    db_path = temp_cache_dir / "test.db"
    db_mgr = DatabaseManager(db_path=db_path)
    db_mgr.init_db()

    assert db_path.exists()


def test_database_caching(temp_cache_dir: Path) -> None:
    """Test setting and getting cached compound."""
    db_path = temp_cache_dir / "test.db"
    db_mgr = DatabaseManager(db_path=db_path)
    db_mgr.init_db()

    compound = PubChemCompound(
        cid=2244,
        input_query="CC(=O)OC1=CC=CC=C1C(=O)O",
        canonical_smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        iupac_name="2-acetoxybenzoic acid",
        molecular_formula="C9H8O4",
        molecular_weight=180.16,
    )

    db_mgr.cache_compound("aspirin", compound)
    cached = db_mgr.get_cached_compound("aspirin")

    assert cached is not None
    assert cached.cid == 2244
    assert cached.iupac_name == "2-acetoxybenzoic acid"