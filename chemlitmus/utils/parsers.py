"""Input parsers for handling various file formats (CSV, TSV, XLSX, SMI, SDF)."""

from pathlib import Path
from typing import List
import pandas as pd

from chemlitmus.logging_config import logger

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


KNOWN_HEADERS = {
    "smiles", "canonical_smiles", "isomeric_smiles", "structure", "compound", "compounds",
    "molecule", "mol", "name", "names", "id", "ids", "cid", "query", "queries", "identifier",
    "inchikey", "inchi_key", "chembl_id", "chebi_id", "kegg_id",
}
_PRIORITY = ["smiles", "canonical_smiles", "isomeric_smiles", "structure", "compound", "molecule",
             "query", "identifier", "inchikey", "inchi_key", "name", "id", "cid"]


def _detect_smiles_column(df: pd.DataFrame) -> str | int:
    """Pick the column holding structures/queries: a recognised header name, else the first column."""
    lowered = {str(c).lower().strip(): c for c in df.columns}
    for name in _PRIORITY:
        if name in lowered:
            return lowered[name]
    return df.columns[0]


def _has_header(first_row: list[str]) -> bool:
    """Decide whether a delimited file's first row is a header.

    A row is treated as a header only when at least one cell is a recognised header name.
    Anything else -- a SMILES string, a compound name, a database identifier -- is data.
    This keeps single-column lists of names or IDs intact (the first entry used to be eaten).
    """
    return any(str(c).lower().strip() in KNOWN_HEADERS for c in first_row)


def _read_delimited(file_path: Path, sep: str) -> List[str]:
    with open(file_path, encoding="utf-8-sig") as fh:
        first = fh.readline()
    cells = [c.strip().strip('"') for c in first.rstrip("\r\n").split(sep)] if first.strip() else []
    if _has_header(cells):
        df = pd.read_csv(file_path, sep=sep, dtype=str, keep_default_na=False)
        col = _detect_smiles_column(df)
    else:
        df = pd.read_csv(file_path, sep=sep, header=None, dtype=str, keep_default_na=False)
        col = 0
    return [s for s in df[col].astype(str).tolist() if s.strip() and s.strip().lower() != "nan"]


def parse_compounds_file(file_path: Path) -> List[str]:
    """
    Parse a file and extract a list of SMILES or CIDs.
    Supports .csv, .tsv, .xlsx, .smi, and .sdf formats.
    """
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    ext = file_path.suffix.lower()
    smiles_list = []

    try:
        if ext in [".csv", ".txt"]:
            smiles_list = _read_delimited(file_path, ",")

        elif ext == ".tsv":
            smiles_list = _read_delimited(file_path, "\t")

        elif ext in [".xlsx", ".xls"]:
            df = pd.read_excel(file_path, dtype=str)
            smiles_col = _detect_smiles_column(df)
            smiles_list = df[smiles_col].dropna().astype(str).tolist()

        elif ext == ".smi":
            with open(file_path) as f:
                # .smi files usually have SMILES as the first space-separated token
                smiles_list = [line.split()[0].strip() for line in f if line.strip()]

        elif ext == ".sdf":
            if not RDKIT_AVAILABLE:
                raise ImportError("RDKit is required to parse SDF files.")
            supplier = Chem.SDMolSupplier(str(file_path))
            for mol in supplier:
                if mol is not None:
                    smiles_list.append(Chem.MolToSmiles(mol))

        else:
            raise ValueError(f"Unsupported file extension: {ext}")

    except Exception as e:
        logger.error(f"Failed to parse {file_path}: {e}")
        raise

    # Clean the list
    return [s.strip() for s in smiles_list if s.strip()]