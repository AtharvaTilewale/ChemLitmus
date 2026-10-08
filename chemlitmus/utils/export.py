"""Export utilities for saving batch results."""

from pathlib import Path
from typing import List
import pandas as pd

from chemlitmus.core.pubchem import PubChemCompound
from chemlitmus.logging_config import logger


def export_results(results: List[PubChemCompound], output_path: Path, fmt: str) -> None:
    """
    Export a list of PubChemCompound results to the specified format.
    
    Args:
        results: List of compound objects.
        output_path: Destination file path.
        fmt: Format ('csv', 'xlsx', 'json').
    """
    if not results:
        logger.warning("No results to export.")
        return

    # Convert Pydantic models to dictionaries
    data = [res.model_dump() for res in results]
    df = pd.DataFrame(data)

    try:
        fmt = fmt.lower().strip()
        if fmt == "csv":
            df.to_csv(output_path, index=False)
        elif fmt == "xlsx":
            df.to_excel(output_path, index=False)
        elif fmt == "json":
            df.to_json(output_path, orient="records", indent=4)
        else:
            raise ValueError(f"Unsupported export format: {fmt}. Use csv, xlsx, or json.")
            
        logger.info(f"Successfully exported {len(results)} records to {output_path}")
        
    except Exception as e:
        logger.error(f"Failed to export results to {output_path}: {e}")
        raise