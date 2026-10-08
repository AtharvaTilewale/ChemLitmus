"""SQLite caching layer for ChemLitmus (Thread-Safe version)."""

import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

from chemlitmus.config import settings
from chemlitmus.core.pubchem import PubChemCompound
from chemlitmus.logging_config import logger


class DatabaseManager:
    """Manages SQLite database initialization and query caching."""

    def __init__(self, db_path: Optional[Path] = None):
        self.db_path = db_path or settings.db_path

    def get_connection(self) -> sqlite3.Connection:
        """Create database connection with multithreading lock timeout."""
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        # timeout=15 forces concurrent threads to wait up to 15s instead of raising "Database is locked"
        conn = sqlite3.connect(self.db_path, timeout=15.0, check_same_thread=False)
        conn.row_factory = sqlite3.Row
        return conn

    def init_db(self) -> None:
        """Initialize SQLite database schema."""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS compound_cache (
                    query_key TEXT PRIMARY KEY,
                    cid INTEGER,
                    smiles TEXT,
                    canonical_smiles TEXT,
                    iupac_name TEXT,
                    molecular_formula TEXT,
                    molecular_weight REAL,
                    inchi TEXT,
                    inchikey TEXT,
                    data_json TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
                """
            )
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_cid ON compound_cache(cid)")
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_smiles ON compound_cache(canonical_smiles)")
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS provider_cache (
                    source TEXT NOT NULL,
                    query_key TEXT NOT NULL,
                    source_id TEXT,
                    inchikey TEXT,
                    record_json TEXT NOT NULL,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    PRIMARY KEY (source, query_key)
                )
                """
            )
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_provider_inchikey ON provider_cache(source, inchikey)")
            conn.commit()
            logger.info(f"Initialized database schema at {self.db_path}")

    def get_cached_compound(self, query_key: str) -> Optional[PubChemCompound]:
        """Fetch cached compound by query key."""
        if not settings.enable_cache:
            return None

        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute(
                    "SELECT data_json FROM compound_cache WHERE query_key = ?",
                    (query_key.lower().strip(),),
                )
                row = cursor.fetchone()
                if row and row["data_json"]:
                    data = json.loads(row["data_json"])
                    compound = PubChemCompound(**data)
                    # Self-healing: if cached compound has a CID but missing canonical_smiles, treat as miss
                    if compound.canonical_smiles is None and compound.cid is not None:
                        return None
                    return compound
        except Exception as e:
            logger.error(f"Cache lookup failed for '{query_key}': {e}")
        return None

    def cache_compound(self, query_key: str, compound: PubChemCompound) -> None:
        """Cache PubChem lookup result."""
        if not settings.enable_cache or not compound:
            return

        try:
            data_dict = compound.model_dump()
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute(
                    """
                    INSERT OR REPLACE INTO compound_cache (
                        query_key, cid, smiles, canonical_smiles, iupac_name,
                        molecular_formula, molecular_weight, inchi, inchikey, data_json, created_at
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        query_key.lower().strip(),
                        compound.cid,
                        compound.isomeric_smiles,  # <-- FIXED: Using isomeric_smiles instead of the missing .smiles attribute
                        compound.canonical_smiles,
                        compound.iupac_name,
                        compound.molecular_formula,
                        compound.molecular_weight,
                        compound.inchi,
                        compound.inchikey,
                        json.dumps(data_dict),
                        datetime.now(timezone.utc).isoformat(),
                    ),
                )
                conn.commit()
        except Exception as e:
            logger.error(f"Failed to write to cache for '{query_key}': {e}")

    # ---- multi-database provider records -------------------------------------------------

    def get_provider_record(self, source: str, query_key: str) -> Optional[dict]:
        """Return the cached ``CompoundRecord`` payload for ``(source, query_key)`` or ``None``."""
        if not settings.enable_cache:
            return None
        try:
            with self.get_connection() as conn:
                row = conn.execute(
                    "SELECT record_json FROM provider_cache WHERE source = ? AND query_key = ?",
                    (source, query_key),
                ).fetchone()
                if row and row["record_json"]:
                    return json.loads(row["record_json"])
        except Exception as e:
            logger.error(f"Provider cache lookup failed for {source}:'{query_key}': {e}")
        return None

    def cache_provider_record(self, source: str, query_keys: List[str], record: dict) -> None:
        """Store one provider record under every key in ``query_keys``."""
        if not settings.enable_cache or not record or not query_keys:
            return
        try:
            payload = json.dumps(record, default=str)
            now = datetime.now(timezone.utc).isoformat()
            with self.get_connection() as conn:
                conn.executemany(
                    """
                    INSERT OR REPLACE INTO provider_cache
                        (source, query_key, source_id, inchikey, record_json, created_at)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    [(source, k, record.get("source_id"), record.get("inchikey"), payload, now) for k in query_keys],
                )
                conn.commit()
        except Exception as e:
            logger.error(f"Failed to write provider cache for {source}: {e}")

    def provider_cache_stats(self) -> Dict[str, int]:
        """Number of cached provider records per source."""
        try:
            with self.get_connection() as conn:
                rows = conn.execute("SELECT source, COUNT(*) AS n FROM provider_cache GROUP BY source").fetchall()
                return {r["source"]: r["n"] for r in rows}
        except Exception:
            return {}
