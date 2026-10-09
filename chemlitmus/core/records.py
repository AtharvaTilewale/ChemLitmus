"""Record-preserving ingestion.

Every row, SDF block or line of an input file becomes one :class:`Record` with a stable internal
id, its position, the source-provided id, the original structure text, every other field, and a
processing status. Nothing is dropped: empty, malformed and unparseable entries are kept with
their positions and issues, so the accounting ``n_total == n_ok + n_empty + n_invalid + n_error``
holds exactly and every downstream count can be reconciled against the input.

Supported formats: CSV, TSV, XLSX/XLS, SMI (SMILES [whitespace name]), SDF, and Parquet (optional,
needs ``pyarrow``: ``pip install 'chemlitmus[parquet]'``). Column roles
(structure, id, endpoint, units, relation, split, date, source) are selected explicitly or
detected from recognised header names; a multi-column table with no recognised structure column
is a schema error unless the caller asks for the first column (``ambiguous="first"``).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from pydantic import BaseModel, Field

from chemlitmus.core.smiles import SmilesParseError, mol_from_smiles, split_smiles_field

try:
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
except ImportError:  # pragma: no cover
    Chem = None

SCHEMA_VERSION = "1"

STATUSES = ["ok", "empty", "invalid", "unsupported", "error"]

STRUCTURE_HEADERS = ["smiles", "canonical_smiles", "isomeric_smiles", "structure", "molecule", "mol", "compound", "query", "smi"]
ID_HEADERS = ["id", "ids", "record_id", "compound_id", "molecule_id", "mol_id", "cid", "chembl_id", "chebi_id", "kegg_id",
              "identifier", "name", "names", "title", "compound_name", "molecule_name"]
ROLE_HEADERS: Dict[str, List[str]] = {
    "structure": STRUCTURE_HEADERS,
    "id": ID_HEADERS,
    "endpoint": ["endpoint", "value", "activity", "label", "y", "target_value", "measurement", "pchembl_value", "standard_value", "ic50", "ki", "ec50"],
    "units": ["units", "unit", "standard_units"],
    "relation": ["relation", "standard_relation", "operator", "qualifier"],
    "split": ["split", "set", "subset", "fold", "partition"],
    "date": ["date", "year", "publication_date", "document_year", "timestamp"],
    "source": ["source", "dataset", "origin", "document_id", "doc_id", "assay_id", "assay"],
    "target": ["target", "target_id", "target_chembl_id", "protein"],
}


class SchemaError(ValueError):
    """The input's columns cannot be mapped to roles without guidance."""


class Issue(BaseModel):
    code: str = Field(description="Stable machine-readable code, e.g. PARSE_WHITESPACE.")
    severity: str = Field(description="info | warning | error")
    message: str
    field: Optional[str] = None
    evidence: Dict[str, Any] = Field(default_factory=dict)


class Record(BaseModel):
    """One input entry with everything needed to trace it back."""

    record_id: str = Field(description="Internal id, unique within the record set (position-based).")
    position: int = Field(description="0-based data row index (tables, SMI) or molecule index (SDF).")
    source_file: str
    source_id: Optional[str] = Field(None, description="Id from the input (id column, SMI name, SDF title). Need not be unique.")
    structure: Optional[str] = Field(None, description="Original structure text exactly as read (SMILES, or MOL block for SDF).")
    structure_format: str = Field("smiles", description="smiles | molblock")
    parsed_smiles: Optional[str] = Field(None, description="RDKit canonical SMILES of the parsed input, no standardisation applied.")
    status: str = Field("ok", description="ok | empty | invalid | unsupported | error")
    issues: List[Issue] = Field(default_factory=list)
    fields: Dict[str, Any] = Field(default_factory=dict, description="Every other column / SDF property, verbatim.")

    @property
    def ok(self) -> bool:
        return self.status == "ok"


class RecordSet(BaseModel):
    source_file: str
    format: str
    schema_version: str = SCHEMA_VERSION
    columns: List[str] = Field(default_factory=list, description="Column names as read (tables) or property names seen (SDF).")
    roles: Dict[str, str] = Field(default_factory=dict, description="role -> column actually used (structure, id, endpoint, ...).")
    n_total: int = 0
    n_ok: int = 0
    n_empty: int = 0
    n_invalid: int = 0
    n_unsupported: int = 0
    n_error: int = 0
    notes: List[str] = Field(default_factory=list)
    records: List[Record] = Field(default_factory=list)

    def ok_records(self) -> List[Record]:
        return [r for r in self.records if r.status == "ok"]

    def reconcile(self) -> bool:
        return self.n_total == self.n_ok + self.n_empty + self.n_invalid + self.n_unsupported + self.n_error == len(self.records)


# ------------------------------------------------------------------------------------ helpers

def _norm(c: Any) -> str:
    return str(c).strip().lower().replace(" ", "_").lstrip("\ufeff")


def detect_roles(columns: Iterable[Any]) -> Dict[str, str]:
    """Map roles to columns by recognised header names (first match in priority order)."""
    cols = list(columns)
    normed = {_norm(c): c for c in cols}
    out: Dict[str, str] = {}
    used = set()
    for role, names in ROLE_HEADERS.items():
        for n in names:
            if n in normed and normed[n] not in used:
                out[role] = str(normed[n]); used.add(normed[n]); break
    return out


def _parse_structure(rec: Record) -> None:
    """Fill ``parsed_smiles``/``status``/``issues`` from ``rec.structure`` (SMILES text)."""
    s = rec.structure
    if s is None or not str(s).strip() or str(s).strip().lower() == "nan":
        rec.status = "empty"
        rec.issues.append(Issue(code="PARSE_EMPTY", severity="warning", message="empty structure field", field="structure"))
        return
    s = str(s)
    try:
        mol = mol_from_smiles(s)
    except SmilesParseError as exc:
        rec.status = "invalid"
        rec.issues.append(Issue(code="PARSE_WHITESPACE", severity="error", field="structure",
                                message="structure field contains whitespace; RDKit would silently keep only the first token",
                                evidence={"text": s, "first_token": s.split()[0] if s.split() else "", "hint": str(exc)}))
        return
    if mol is None:
        rec.status = "invalid"
        rec.issues.append(Issue(code="PARSE_INVALID", severity="error", field="structure", message="RDKit cannot parse the structure", evidence={"text": s}))
        return
    rec.parsed_smiles = Chem.MolToSmiles(mol)
    rec.status = "ok"


def _finish(rs: RecordSet) -> RecordSet:
    rs.n_total = len(rs.records)
    for r in rs.records:
        setattr(rs, f"n_{r.status}", getattr(rs, f"n_{r.status}") + 1)
    assert rs.reconcile()
    return rs


# ------------------------------------------------------------------------------------ adapters

def _read_table(path: Path, sep: Optional[str], structure_column: Optional[str], id_column: Optional[str],
                roles: Optional[Dict[str, str]], ambiguous: str, sheet: Optional[str]) -> RecordSet:
    import pandas as pd

    ext = path.suffix.lower()
    fmt = {".csv": "csv", ".tsv": "tsv", ".txt": "csv"}.get(ext, "xlsx")
    if ext in (".xlsx", ".xls"):
        df = pd.read_excel(path, dtype=str, keep_default_na=False, sheet_name=sheet or 0)
        has_header = True
    else:
        with open(path, encoding="utf-8-sig") as fh:
            first = fh.readline()
        cells = [c.strip().strip('"') for c in first.rstrip("\r\n").split(sep or ",")] if first.strip() else []
        known = {n for names in ROLE_HEADERS.values() for n in names}
        explicit = {c for c in [structure_column, id_column, *(roles or {}).values()] if c}
        has_header = any(_norm(c) in known for c in cells)
        if not has_header and explicit and explicit <= set(cells):
            has_header = True                       # the caller named these columns: the first row is the header
        if not has_header and len(cells) > 1:
            # A multi-column first row in which no cell is a parseable SMILES is taken as a header
            # (``foo,bar``). A row that does contain a SMILES is data. Single-column files are not
            # guessed at: a lone unrecognised word is a query, not a header.
            def _is_smiles(c: str) -> bool:
                try:
                    return bool(c) and mol_from_smiles(c) is not None
                except SmilesParseError:
                    return False
            if not any(_is_smiles(c) for c in cells):
                has_header = True
        df = pd.read_csv(path, sep=sep or ",", dtype=str, keep_default_na=False, header=0 if has_header else None, encoding="utf-8-sig")
        if not has_header:
            df.columns = [f"col{i}" for i in range(df.shape[1])]
    columns = [str(c) for c in df.columns]
    detected = detect_roles(columns) if has_header else {}
    used_roles = dict(detected)
    if roles:
        for role, col in roles.items():
            if col not in columns:
                raise SchemaError(f"Column {col!r} for role {role!r} not found. Columns: {columns}")
            used_roles[role] = col
    if structure_column:
        if structure_column not in columns:
            raise SchemaError(f"Structure column {structure_column!r} not found. Columns: {columns}")
        used_roles["structure"] = structure_column
    if id_column:
        if id_column not in columns:
            raise SchemaError(f"Id column {id_column!r} not found. Columns: {columns}")
        used_roles["id"] = id_column
    notes: List[str] = []
    if "structure" not in used_roles:
        if len(columns) == 1 or ambiguous == "first":
            used_roles["structure"] = columns[0]
            if len(columns) > 1:
                notes.append(f"no recognised structure column; using the first column {columns[0]!r} because ambiguous='first'")
        else:
            raise SchemaError(
                f"Cannot tell which of {columns} holds the structures. Pass structure_column=... "
                f"(recognised header names: {', '.join(STRUCTURE_HEADERS)})."
            )
    if not has_header:
        notes.append("no header row recognised; columns named col0, col1, ...")
    rs = RecordSet(source_file=str(path), format=fmt, columns=columns, roles=used_roles, notes=notes)
    scol, idcol = used_roles["structure"], used_roles.get("id")
    for i, row in enumerate(df.itertuples(index=False, name=None)):
        d = dict(zip(columns, row))
        rec = Record(record_id=f"r{i:06d}", position=i, source_file=str(path),
                     source_id=(str(d[idcol]).strip() or None) if idcol else None,
                     structure=(None if str(d.get(scol, "")).strip() == "" else str(d[scol])),
                     fields={k: v for k, v in d.items() if k != scol})
        _parse_structure(rec)
        rs.records.append(rec)
    return _finish(rs)


def _read_smi(path: Path) -> RecordSet:
    rs = RecordSet(source_file=str(path), format="smi", columns=["smiles", "name"], roles={"structure": "smiles", "id": "name"})
    with open(path, encoding="utf-8-sig") as fh:
        lines = fh.read().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("#"):
            rec = Record(record_id=f"r{i:06d}", position=i, source_file=str(path), structure=None, status="unsupported",
                         issues=[Issue(code="SMI_COMMENT", severity="info", message="comment line")], fields={"raw": line})
            rs.records.append(rec); continue
        smi, name = split_smiles_field(line)
        rec = Record(record_id=f"r{i:06d}", position=i, source_file=str(path), source_id=name, structure=smi or None,
                     fields={"name": name} if name else {})
        _parse_structure(rec)
        rs.records.append(rec)
    return _finish(rs)


_PROP_RE = re.compile(r"^>\s*<([^>]+)>[^\n]*\n(.*?)(?=\n>\s*<|\n\$\$\$\$|\Z)", re.S | re.M)


def _read_sdf(path: Path) -> RecordSet:
    rs = RecordSet(source_file=str(path), format="sdf", roles={"structure": "molblock", "id": "title"})
    supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=True)
    props_seen: List[str] = []
    for i in range(len(supplier)):
        raw = supplier.GetItemText(i)
        title = raw.split("\n", 1)[0].strip() if raw else ""
        props = {m.group(1).strip(): m.group(2).strip() for m in _PROP_RE.finditer(raw or "")}
        for k in props:
            if k not in props_seen:
                props_seen.append(k)
        rec = Record(record_id=f"r{i:06d}", position=i, source_file=str(path), source_id=title or None,
                     structure=raw, structure_format="molblock", fields=props)
        mol = supplier[i]
        if mol is None:
            rec.status = "invalid"
            rec.issues.append(Issue(code="SDF_PARSE_FAILED", severity="error", field="structure",
                                    message="RDKit could not read or sanitise this SDF record", evidence={"title": title}))
        else:
            rec.parsed_smiles = Chem.MolToSmiles(mol)
            rec.status = "ok"
        rs.records.append(rec)
    rs.columns = ["title"] + props_seen
    return _finish(rs)


def _read_parquet(path: Path, structure_column, id_column, roles, ambiguous) -> RecordSet:
    """Optional format: needs ``pyarrow`` (``pip install chemlitmus[parquet]``)."""
    import pandas as pd

    try:
        df = pd.read_parquet(path)
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError("Reading Parquet needs pyarrow: pip install 'chemlitmus[parquet]'") from exc
    df = df.astype(str).fillna("")
    columns = [str(c) for c in df.columns]
    used = dict(detect_roles(columns))
    for role, col in (roles or {}).items():
        if col not in columns:
            raise SchemaError(f"Column {col!r} for role {role!r} not found. Columns: {columns}")
        used[role] = col
    for role, col in (("structure", structure_column), ("id", id_column)):
        if col:
            if col not in columns:
                raise SchemaError(f"Column {col!r} not found. Columns: {columns}")
            used[role] = col
    notes: List[str] = []
    if "structure" not in used:
        if len(columns) == 1 or ambiguous == "first":
            used["structure"] = columns[0]
            if len(columns) > 1:
                notes.append(f"no recognised structure column; using {columns[0]!r} because ambiguous='first'")
        else:
            raise SchemaError(f"Cannot tell which of {columns} holds the structures. Pass structure_column=...")
    rs = RecordSet(source_file=str(path), format="parquet", columns=columns, roles=used, notes=notes)
    scol, idcol = used["structure"], used.get("id")
    for i, row in enumerate(df.itertuples(index=False, name=None)):
        d = dict(zip(columns, row))
        rec = Record(record_id=f"r{i:06d}", position=i, source_file=str(path),
                     source_id=(str(d[idcol]).strip() or None) if idcol else None,
                     structure=(None if str(d.get(scol, "")).strip() == "" else str(d[scol])),
                     fields={k: v for k, v in d.items() if k != scol})
        _parse_structure(rec)
        rs.records.append(rec)
    return _finish(rs)


def read_records(
    path: str | Path,
    structure_column: Optional[str] = None,
    id_column: Optional[str] = None,
    roles: Optional[Dict[str, str]] = None,
    ambiguous: str = "error",
    sheet: Optional[str] = None,
) -> RecordSet:
    """Read any supported file into a :class:`RecordSet`.

    Args:
        path: CSV/TSV/TXT/XLSX/XLS/SMI/SDF/Parquet file.
        structure_column: Column holding SMILES (tables). Overrides header detection.
        id_column: Column holding the source id.
        roles: Extra role -> column assignments (endpoint, units, relation, split, date, source, target).
        ambiguous: ``"error"`` (default) raises :class:`SchemaError` when a multi-column table has no
            recognised structure column; ``"first"`` uses the first column and records a note.
        sheet: Worksheet name for Excel files.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {p}")
    ext = p.suffix.lower()
    if ext in (".csv", ".txt", ".xlsx", ".xls"):
        return _read_table(p, None, structure_column, id_column, roles, ambiguous, sheet)
    if ext == ".tsv":
        return _read_table(p, "\t", structure_column, id_column, roles, ambiguous, sheet)
    if ext == ".smi":
        return _read_smi(p)
    if ext == ".sdf":
        return _read_sdf(p)
    if ext in (".parquet", ".pq"):
        return _read_parquet(p, structure_column, id_column, roles, ambiguous)
    raise ValueError(f"Unsupported file extension: {ext}")


__all__ = ["Record", "RecordSet", "Issue", "SchemaError", "read_records", "detect_roles", "STATUSES", "ROLE_HEADERS", "SCHEMA_VERSION"]
