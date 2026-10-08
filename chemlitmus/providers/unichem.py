"""UniChem cross-references: one InChIKey → identifiers in ~40 databases."""

from __future__ import annotations

import threading
import time
from typing import Dict

import requests

from chemlitmus.providers.base import ProviderError

_URL = "https://www.ebi.ac.uk/unichem/rest/inchikey/{}"

# UniChem src_id -> short name. Only the sources a chemist is likely to want are named; the rest
# are reported as "unichem_src_<id>" so nothing is dropped silently.
UNICHEM_SOURCES: Dict[str, str] = {
    "1": "chembl", "2": "drugbank", "3": "pdbe", "4": "gtopdb", "5": "pubchem_dotf", "6": "kegg_ligand",
    "7": "chebi", "9": "zinc", "10": "emolecules", "14": "fdasrs", "15": "surechembl", "17": "pharmgkb",
    "18": "hmdb", "20": "selleck", "21": "pubchem_tpharma", "22": "pubchem", "23": "mcule", "24": "nmrshiftdb2",
    "25": "lincs", "26": "actor", "27": "recon", "28": "molport", "29": "nikkaji", "31": "bindingdb",
    "32": "comptox", "33": "lipidmaps", "34": "drugcentral", "35": "carotenoiddb", "36": "metabolights",
    "37": "brenda", "38": "rhea", "39": "chemicalbook", "41": "swisslipids", "45": "dailymed", "46": "clinicaltrials",
    "47": "rxnorm", "48": "medchemexpress", "49": "probes_drugs", "50": "mlsmr", "51": "foodb", "84": "cas",
}

_lock = threading.Lock()
_last = 0.0


def unichem_xrefs(inchikey: str, timeout: float = 15.0) -> Dict[str, str]:
    """Return ``{source_name: identifier}`` for every UniChem source holding this InChIKey."""
    global _last
    with _lock:
        wait = 0.34 - (time.time() - _last)
        if wait > 0:
            time.sleep(wait)
        _last = time.time()
    try:
        r = requests.get(_URL.format(inchikey), timeout=timeout,
                         headers={"User-Agent": "chemlitmus (https://github.com/AtharvaTilewale/ChemLitmus)"})
    except requests.RequestException as exc:
        raise ProviderError(f"UniChem: network error: {exc}") from exc
    if r.status_code == 404:
        return {}
    if r.status_code != 200:
        raise ProviderError(f"UniChem: HTTP {r.status_code}")
    try:
        rows = r.json()
    except ValueError as exc:
        raise ProviderError("UniChem: non-JSON response") from exc
    if isinstance(rows, dict) and rows.get("error"):
        return {}
    out: Dict[str, str] = {}
    for row in rows if isinstance(rows, list) else []:
        sid, cid = str(row.get("src_id", "")), str(row.get("src_compound_id", ""))
        if not sid or not cid:
            continue
        name = UNICHEM_SOURCES.get(sid, f"unichem_src_{sid}")
        out.setdefault(name, cid)
    return out


__all__ = ["unichem_xrefs", "UNICHEM_SOURCES"]
