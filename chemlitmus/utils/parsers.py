"""Input parsers (compatibility layer over :mod:`chemlitmus.core.records`)."""

from pathlib import Path
from typing import List

def parse_compounds_file(file_path: Path) -> List[str]:
    """Return the structure strings of a file as a flat list (compatibility wrapper).

    This is the historical interface: a list of SMILES (or query strings) with blanks removed and,
    for SDF, the canonical SMILES of the records RDKit could read. It **discards** ids, other
    columns and failed records. New code should use :func:`chemlitmus.core.records.read_records`,
    which keeps every record with its position, metadata and status.
    """
    from chemlitmus.core.records import read_records

    rs = read_records(file_path, ambiguous="first")
    out: List[str] = []
    for r in rs.records:
        if r.status == "empty" or r.status == "unsupported":
            continue
        if r.structure_format == "molblock":
            if r.parsed_smiles:
                out.append(r.parsed_smiles)
        elif r.structure is not None:
            out.append(str(r.structure).strip())
    return out
