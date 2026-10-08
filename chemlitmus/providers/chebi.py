"""ChEBI provider (EBI ChEBI backend API)."""

from __future__ import annotations

import re
from typing import Any, Dict, Optional

from chemlitmus.providers.base import CompoundRecord, Provider, _float, _int

_BASE = "https://www.ebi.ac.uk/chebi/backend/api/public"
_ID_RE = re.compile(r"^(CHEBI:)?\d+$", re.I)


_TAG = re.compile(r"<[^>]+>")


def _untag(s):
    """ChEBI names carry presentational HTML (``<small>L</small>-thyroxine``, ``<sup>``)."""
    return _TAG.sub("", s) if isinstance(s, str) else s


class ChEBIProvider(Provider):
    key = "chebi"
    name = "ChEBI"
    homepage = "https://www.ebi.ac.uk/chebi"
    min_interval = 0.34

    def looks_like_id(self, query: str) -> bool:
        return query.upper().startswith("CHEBI:")

    @staticmethod
    def _numeric(identifier: str) -> str:
        return identifier.split(":")[-1]

    def _to_record(self, c: Dict[str, Any], query: str) -> CompoundRecord:
        acc = c.get("chebi_accession") or f"CHEBI:{c.get('id')}"
        chem = c.get("chemical_data") or {}
        st = c.get("default_structure") or {}
        names = c.get("names") or {}
        synonyms = []
        for group in ("SYNONYM", "IUPAC NAME", "INN", "BRAND NAME"):
            for n in names.get(group, []) if isinstance(names, dict) else []:
                v = n.get("name") if isinstance(n, dict) else n
                if v:
                    synonyms.append(v)
        xrefs: Dict[str, str] = {}
        for group, entries in (c.get("database_accessions") or {}).items():
            for e in entries:
                acc_no = str(e.get("accession_number") or "")
                etype = (e.get("type") or group or "").upper()
                src = (e.get("source_name") or "").lower()
                if not acc_no:
                    continue
                if etype == "CAS":
                    xrefs.setdefault("cas", acc_no)
                elif etype == "MANUAL_X_REF" and src in ("drugbank", "pdbechem", "drugcentral", "hmdb", "lipid maps", "kegg compound", "pubchem"):
                    key = {"pdbechem": "pdbe", "lipid maps": "lipidmaps", "kegg compound": "kegg"}.get(src, src)
                    xrefs.setdefault(key, acc_no)
        return CompoundRecord(
            source=self.key, source_id=acc, query=query,
            url=f"{self.homepage}/searchId.do?chebiId={acc}",
            name=_untag(c.get("name") or c.get("ascii_name")), synonyms=[_untag(s) for s in synonyms[:25]],
            smiles=st.get("smiles"), inchi=st.get("standard_inchi"), inchikey=st.get("standard_inchi_key"),
            formula=chem.get("formula"),
            molecular_weight=_float(chem.get("mass")), monoisotopic_mass=_float(chem.get("monoisotopic_mass")),
            charge=_int(chem.get("charge")),
            extra={k: v for k, v in {"stars": c.get("stars"), "definition": c.get("definition")}.items() if v is not None},
            cross_refs=xrefs,
        )

    def by_id(self, identifier: str) -> Optional[CompoundRecord]:
        c = self._get(f"{_BASE}/compound/{self._numeric(identifier)}/")
        return self._to_record(c, identifier) if c and c.get("chebi_accession") else None

    def _search(self, term: str, query: str, exact_name: bool = False) -> Optional[CompoundRecord]:
        data = self._get(f"{_BASE}/es_search/", {"term": term, "size": 10})
        hits = [h.get("_source", {}) for h in ((data or {}).get("results") or [])]
        if not hits:
            return None
        lterm = term.strip().lower()
        if exact_name:
            # 1) a hit whose primary name is the term; 2) whose synonym list contains it
            for h in hits:
                if _untag(h.get("name") or "").lower() == lterm or (h.get("ascii_name") or "").lower() == lterm:
                    return self.by_id(h["chebi_accession"])
            for h in sorted(hits, key=lambda x: -(x.get("stars") or 0)):
                rec = self.by_id(h["chebi_accession"])
                if rec and any(s.lower() == lterm for s in rec.synonyms):
                    return rec
            # 3) fall back to the highest-ranked curated entry
        hits.sort(key=lambda h: (-(h.get("stars") or 0),))
        acc = hits[0].get("chebi_accession")
        return self.by_id(acc) if acc else None

    def by_inchikey(self, inchikey: str) -> Optional[CompoundRecord]:
        data = self._get(f"{_BASE}/es_search/", {"term": inchikey, "size": 10})
        for h in ((data or {}).get("results") or []):
            src = h.get("_source", {})
            if (src.get("inchikey") or "").upper() == inchikey.upper():
                return self.by_id(src["chebi_accession"])
        return None

    def by_name(self, name: str) -> Optional[CompoundRecord]:
        return self._search(name, name, exact_name=True)
