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


def _detect_smiles_column(df: pd.DataFrame) -> str | int:
    """Auto-detect the SMILES column in a DataFrame."""
    common_names = ["smiles", "canonical_smiles", "structure", "compound"]
    
    # Check if any column header matches common names
    for col in df.columns:
        if str(col).lower().strip() in common_names:
            return col
            
    # Fallback: Assume it's a headless CSV and the first column is SMILES
    return df.columns[0]


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
            # Read first row to guess if header exists
            df_test = pd.read_csv(file_path, nrows=1)
            smiles_col = _detect_smiles_column(df_test)
            
            # If the detected column is literally a SMILES string (e.g., 'c1ccccc1'), it has no header
            if any(c in str(smiles_col) for c in ["=", "#", "(", ")", "c", "C"]):
                df = pd.read_csv(file_path, header=None)
                smiles_list = df[0].dropna().astype(str).tolist()
            else:
                df = pd.read_csv(file_path)
                smiles_list = df[smiles_col].dropna().astype(str).tolist()

        elif ext == ".tsv":
            df = pd.read_csv(file_path, sep="\t")
            smiles_col = _detect_smiles_column(df)
            smiles_list = df[smiles_col].dropna().astype(str).tolist()

        elif ext in [".xlsx", ".xls"]:
            df = pd.read_excel(file_path)
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