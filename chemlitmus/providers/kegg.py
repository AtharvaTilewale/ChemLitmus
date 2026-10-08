"""KEGG COMPOUND provider (KEGG REST, flat-file format)."""

from __future__ import annotations

import re
from typing import Dict, List, Optional

from chemlitmus.providers.base import CompoundRecord, Provider, _float

_BASE = "https://rest.kegg.jp"
_ID_RE = re.compile(r"^(cpd:)?C\d{5}$", re.I)


def _parse_flat(text: str) -> Dict[str, List[str]]:
    """Parse KEGG's 12-column flat-file format into {FIELD: [lines]}."""
    out: Dict[str, List[str]] = {}
    field = None
    for line in text.splitlines():
        if line.startswith("///"):
            break
        if line[:12].strip():
            field = line[:12].strip()
            out.setdefault(field, [])
        if field is not None:
            val = line[12:].strip()
            if val:
                out[field].append(val)
    return out


class KEGGProvider(Provider):
    key = "kegg"
    name = "KEGG"
    homepage = "https://www.kegg.jp"
    min_interval = 0.4                                     # KEGG asks for ≤3 requests/s

    def looks_like_id(self, query: str) -> bool:
        return bool(_ID_RE.match(query))

    def by_id(self, identifier: str) -> Optional[CompoundRecord]:
        cid = identifier.split(":")[-1].upper()
        text = self._get(f"{_BASE}/get/{cid}", json=False)
        if not text:
            return None
        f = _parse_flat(text)
        names = [n.rstrip(";") for n in f.get("NAME", [])]
        xrefs: Dict[str, str] = {}
        for ln in f.get("DBLINKS", []):
            if ":" in ln:
                src, val = ln.split(":", 1)
                src = src.strip().lower(); val = val.strip().split()[0]
                if src in ("pubchem", "chebi", "chembl", "cas", "lipidmaps", "hmdb", "drugbank", "knapsack"):
                    xrefs[src] = f"CHEBI:{val}" if src == "chebi" and not val.upper().startswith("CHEBI") else val
        # KEGG's PubChem link is a *substance* (SID) not a compound (CID); label it so
        if "pubchem" in xrefs:
            xrefs["pubchem_sid"] = xrefs.pop("pubchem")
        return CompoundRecord(
            source=self.key, source_id=cid, query=identifier,
            url=f"{self.homepage}/entry/{cid}",
            name=names[0] if names else None, synonyms=names[1:25],
            formula=(f.get("FORMULA") or [None])[0],
            molecular_weight=_float((f.get("MOL_WEIGHT") or [None])[0]),
            monoisotopic_mass=_float((f.get("EXACT_MASS") or [None])[0]),
            extra={k: v for k, v in {
                "pathways": [p.split()[0] for p in f.get("PATHWAY", [])][:20],
                "enzymes": " ".join(f.get("ENZYME", [])).split()[:20] or None,
                "brite": f.get("BRITE", [None])[0],
            }.items() if v},
            cross_refs=xrefs,
        )

    def by_name(self, name: str) -> Optional[CompoundRecord]:
        text = self._get(f"{_BASE}/find/compound/{name}", json=False)
        if not text or not text.strip():
            return None
        lname = name.lower()
        best = None
        for line in text.splitlines():
            cid, _, names = line.partition("\t")
            cands = [n.strip().lower() for n in names.split(";")]
            if lname in cands:
                best = cid; break
            best = best or cid
        return self.by_id(best) if best else None

    def by_inchikey(self, inchikey: str) -> Optional[CompoundRecord]:
        # KEGG has no InChIKey endpoint; resolve through ChEBI's cross-references when available.
        try:
            from chemlitmus.providers.chebi import ChEBIProvider
            ch = ChEBIProvider().by_inchikey(inchikey)
        except Exception:
            ch = None
        kid = (ch.cross_refs.get("kegg") if ch else None)
        return self.by_id(kid) if kid else None
