"""Common schema and base class for compound-database providers.

Every provider turns a query (name, SMILES, InChIKey, or a source-specific identifier) into a
:class:`CompoundRecord` — the same shape regardless of where it came from — so results from
PubChem, ChEMBL, ChEBI and KEGG can be compared, merged and tabulated together.
"""

from __future__ import annotations

import threading
import time
from datetime import datetime, timezone
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional

import requests
from pydantic import BaseModel, Field
from requests.adapters import HTTPAdapter

from chemlitmus.logging_config import logger


QUERY_TYPES: List[str] = ["auto", "name", "smiles", "inchikey", "id"]
"""Query types every provider understands. ``id`` is the provider's own identifier
(PubChem CID, ChEMBL ID, ChEBI ID, KEGG compound ID)."""


class CompoundRecord(BaseModel):
    """One compound as described by one database."""

    source: str = Field(description="Provider key, e.g. 'pubchem', 'chembl', 'chebi', 'kegg'.")
    source_id: str = Field(description="The provider's identifier for this record.")
    url: Optional[str] = Field(None, description="Landing page for the record.")
    query: str = Field(description="The query string that produced this record.")

    name: Optional[str] = Field(None, description="Preferred / title name.")
    synonyms: List[str] = Field(default_factory=list)
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    monoisotopic_mass: Optional[float] = None
    charge: Optional[int] = None

    # Optional computed properties where the source supplies them
    xlogp: Optional[float] = None
    hbd: Optional[int] = None
    hba: Optional[int] = None
    tpsa: Optional[float] = None
    rotatable_bonds: Optional[int] = None

    # Source-specific annotations that are worth carrying but do not fit the common schema
    extra: Dict[str, Any] = Field(default_factory=dict, description="Provider-specific fields (ChEMBL max_phase, ChEBI stars, KEGG pathways …).")
    cross_refs: Dict[str, str] = Field(default_factory=dict, description="Identifiers in other databases, keyed by provider.")

    # Retrieval provenance
    retrieved_at: Optional[str] = Field(None, description="UTC ISO-8601 timestamp of the network fetch that produced this record.")
    from_cache: bool = Field(False, description="True when this copy came from the local SQLite cache rather than the network.")
    cached_at: Optional[str] = Field(None, description="When the cached copy was first retrieved; equal to retrieved_at on a fresh fetch. Lets you tell stale data from fresh.")
    source_version: Optional[str] = Field(None, description="Release or version string of the source record, where the database supplies one.")
    field_provenance: Dict[str, str] = Field(default_factory=dict, description="Merged records only: field name -> the source it was taken from.")

    def summary_line(self) -> str:
        bits = [f"{self.source}:{self.source_id}"]
        if self.name:
            bits.append(self.name)
        if self.formula:
            bits.append(self.formula)
        return "  ".join(bits)


class ProviderError(RuntimeError):
    """A provider could not answer (network, HTTP error, unexpected payload)."""


class Provider(ABC):
    """Base class: a rate-limited HTTP session plus the query-dispatch contract."""

    key: str = ""
    name: str = ""
    homepage: str = ""
    min_interval: float = 0.34          # seconds between requests (≈3 req/s, polite default)
    timeout: float = 15.0
    retries: int = 2
    backoff: float = 0.5

    _locks: Dict[str, threading.Lock] = {}
    _last: Dict[str, float] = {}
    _local = threading.local()

    def __init__(self) -> None:
        self.session = requests.Session()
        adapter = HTTPAdapter(max_retries=0)       # retries are handled in _get, under the deadline
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)
        self.session.headers["User-Agent"] = "chemlitmus (https://github.com/AtharvaTilewale/ChemLitmus)"
        Provider._locks.setdefault(self.key, threading.Lock())
        Provider._last.setdefault(self.key, 0.0)

    # ---- plumbing ------------------------------------------------------------------------

    def _throttle(self) -> None:
        with Provider._locks[self.key]:
            wait = self.min_interval - (time.time() - Provider._last[self.key])
            if wait > 0:
                time.sleep(wait)
            Provider._last[self.key] = time.time()

    def _get(self, url: str, params: Optional[Dict[str, Any]] = None, *, json: bool = True) -> Any:
        """GET with rate limiting, deadline-aware retries on transient errors, and 404 -> None."""
        deadline = getattr(self._local, "deadline", None)
        last_err: Optional[str] = None
        for attempt in range(self.retries + 1):
            if deadline is not None and time.time() >= deadline:
                raise ProviderError(f"{self.name}: lookup deadline exceeded" + (f" ({last_err})" if last_err else ""))
            timeout = self.timeout if deadline is None else max(1.0, min(self.timeout, deadline - time.time()))
            self._throttle()
            logger.debug("%s GET %s %s (attempt %d)", self.key, url, params or "", attempt + 1)
            try:
                r = self.session.get(url, params=params, timeout=timeout)
            except requests.RequestException as exc:
                last_err = f"network error: {exc.__class__.__name__}"
                r = None
            if r is not None:
                if r.status_code == 404:
                    return None
                if r.status_code == 200:
                    if not json:
                        return r.text
                    try:
                        return r.json()
                    except ValueError as exc:
                        raise ProviderError(f"{self.name}: non-JSON response from {url}") from exc
                if r.status_code not in (429, 500, 502, 503, 504):
                    raise ProviderError(f"{self.name}: HTTP {r.status_code}")
                detail = ""
                try:
                    fault = r.json().get("Fault", {})
                    detail = f" ({fault.get('Code')}: {fault.get('Message')})" if fault else ""
                except ValueError:
                    pass
                last_err = f"HTTP {r.status_code}{detail}"
            if attempt < self.retries:
                sleep = self.backoff * (2 ** attempt)
                if deadline is not None:
                    sleep = min(sleep, max(0.0, deadline - time.time()))
                time.sleep(sleep)
        raise ProviderError(f"{self.name}: {last_err or 'request failed'}")

    # ---- contract ------------------------------------------------------------------------

    @abstractmethod
    def by_inchikey(self, inchikey: str) -> Optional[CompoundRecord]: ...

    @abstractmethod
    def by_name(self, name: str) -> Optional[CompoundRecord]: ...

    @abstractmethod
    def by_id(self, identifier: str) -> Optional[CompoundRecord]: ...

    def by_smiles(self, smiles: str) -> Optional[CompoundRecord]:
        """Default: canonicalise to an InChIKey locally and look that up. Providers with a
        native structure search override this."""
        ik = _smiles_to_inchikey(smiles)
        return self.by_inchikey(ik) if ik else None

    @abstractmethod
    def looks_like_id(self, query: str) -> bool:
        """Whether ``query`` has the shape of this provider's native identifier."""

    def lookup(
        self,
        query: str,
        query_type: str = "auto",
        deadline: Optional[float] = None,
        use_cache: bool = True,
    ) -> Optional[CompoundRecord]:
        """Dispatch by query type, consulting the local SQLite cache first.

        ``deadline`` is an absolute ``time.time()`` after which no further HTTP requests are
        started and in-flight ones are capped. A found record is cached under the query, its
        native identifier and its InChIKey, so later lookups by any of them are offline.
        """
        qt = query_type.lower()
        if qt not in QUERY_TYPES:
            raise ValueError(f"Unknown query type {query_type!r}. Valid: {QUERY_TYPES}")
        key = self._cache_key(query, qt)
        if use_cache:
            cached = _cache_get(self.key, key)
            if cached is not None:
                cached.from_cache = True
                return cached
        self._local.deadline = deadline
        try:
            rec = self._lookup(query, query_type)
        finally:
            self._local.deadline = None
        if rec is not None:
            now = datetime.now(timezone.utc).isoformat(timespec="seconds")
            rec.retrieved_at = now
            rec.cached_at = rec.cached_at or now
            rec.from_cache = False
            if use_cache:
                _cache_put(self.key, rec, extra_keys=[key])
        return rec

    def _cache_key(self, query: str, query_type: str) -> str:
        """Normalised cache key: ``id:``/``inchikey:``/``smiles:``/``text:`` prefix plus value."""
        q = query.strip()
        if query_type == "auto":
            if _looks_like_inchikey(q):
                query_type = "inchikey"
            elif self.looks_like_id(q):
                query_type = "id"
            elif _is_valid_smiles(q):
                query_type = "smiles"
            else:
                query_type = "text"
        elif query_type == "name":
            query_type = "text"
        if query_type == "smiles":
            return f"smiles:{_canonical_smiles(q) or q}"
        return f"{query_type}:{q.lower()}"

    def _lookup(self, query: str, query_type: str) -> Optional[CompoundRecord]:
        q = query.strip()
        qt = query_type.lower()
        if qt not in QUERY_TYPES:
            raise ValueError(f"Unknown query type {query_type!r}. Valid: {QUERY_TYPES}")
        if qt == "inchikey" or (qt == "auto" and _looks_like_inchikey(q)):
            return self.by_inchikey(q)
        if qt == "id" or (qt == "auto" and self.looks_like_id(q)):
            return self.by_id(q)
        if qt == "smiles" or (qt == "auto" and _is_valid_smiles(q)):
            return self.by_smiles(q)
        return self.by_name(q)


# ---- shared helpers ----------------------------------------------------------------------

_DB = None
_DB_LOCK = threading.Lock()


def _db():
    """Lazily opened shared cache handle (same SQLite file as the PubChem cache).

    Providers run in parallel threads, so creation and schema initialisation are serialised.
    """
    global _DB
    with _DB_LOCK:
        if _DB is None:
            from chemlitmus.core.database import DatabaseManager
            db = DatabaseManager()
            db.init_db()
            _DB = db
    return _DB


def _cache_get(source: str, key: str) -> Optional[CompoundRecord]:
    from chemlitmus.config import settings
    if not settings.enable_cache:
        return None
    data = _db().get_provider_record(source, key)
    if data is None:
        return None
    try:
        return CompoundRecord.model_validate(data)
    except Exception:  # corrupt or outdated payload: treat as a miss
        return None


def _cache_put(source: str, rec: CompoundRecord, extra_keys: Optional[List[str]] = None) -> None:
    from chemlitmus.config import settings
    if not settings.enable_cache:
        return
    keys = set(extra_keys or [])
    if rec.source_id:
        keys.add(f"id:{rec.source_id.lower()}")
    if rec.inchikey:
        keys.add(f"inchikey:{rec.inchikey.lower()}")
    _db().cache_provider_record(source, sorted(keys), rec.model_dump(mode="json"))


def _canonical_smiles(s: str) -> Optional[str]:
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(s)
        return Chem.MolToSmiles(m) if m is not None else None
    except ImportError:  # pragma: no cover
        return None


def _looks_like_inchikey(s: str) -> bool:
    return len(s) == 27 and s[14] == "-" and s[25] == "-" and s.replace("-", "").isalnum() and s.isupper()


def _is_valid_smiles(s: str) -> bool:
    if any(ch.isspace() for ch in s):
        return False
    try:
        from rdkit import Chem
        return Chem.MolFromSmiles(s) is not None
    except ImportError:  # pragma: no cover
        return False


def _smiles_to_inchikey(s: str) -> Optional[str]:
    try:
        from rdkit import Chem
        from rdkit.Chem import inchi
        m = Chem.MolFromSmiles(s)
        return inchi.MolToInchiKey(m) if m is not None else None
    except ImportError:  # pragma: no cover
        return None


def _float(v: Any) -> Optional[float]:
    try:
        return float(v) if v not in (None, "") else None
    except (TypeError, ValueError):
        return None


def _int(v: Any) -> Optional[int]:
    try:
        return int(float(v)) if v not in (None, "") else None
    except (TypeError, ValueError):
        return None


__all__ = ["QUERY_TYPES", "CompoundRecord", "Provider", "ProviderError"]
