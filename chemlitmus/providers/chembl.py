"""ChEMBL provider (EBI web services)."""

from __future__ import annotations

import re
from typing import Any, Dict, Optional

from chemlitmus.providers.base import CompoundRecord, Provider, _float, _int

_BASE = "https://www.ebi.ac.uk/chembl/api/data"
_ID_RE = re.compile(r"^CHEMBL\d+$", re.I)


class ChEMBLProvider(Provider):
    key = "chembl"
    name = "ChEMBL"
    homepage = "https://www.ebi.ac.uk/chembl"
    min_interval = 0.25

    def looks_like_id(self, query: str) -> bool:
        return bool(_ID_RE.match(query))

    def _to_record(self, m: Dict[str, Any], query: str) -> CompoundRecord:
        st = m.get("molecule_structures") or {}
        pr = m.get("molecule_properties") or {}
        cid = m.get("molecule_chembl_id", "")
        synonyms = sorted({s.get("molecule_synonym") for s in (m.get("molecule_synonyms") or []) if s.get("molecule_synonym")})
        extra = {k: m.get(k) for k in ("max_phase", "molecule_type", "first_approval", "oral", "parenteral", "topical",
                                        "withdrawn_flag", "atc_classifications", "natural_product", "prodrug", "chirality")}
        extra["qed_weighted"] = _float(pr.get("qed_weighted"))
        extra["num_ro5_violations"] = _int(pr.get("num_ro5_violations"))
        extra["molecular_species"] = pr.get("molecular_species")
        return CompoundRecord(
            source=self.key, source_id=cid, query=query,
            url=f"{self.homepage}/compound_report_card/{cid}/" if cid else None,
            name=m.get("pref_name"), synonyms=synonyms[:25],
            smiles=st.get("canonical_smiles"), inchi=st.get("standard_inchi"), inchikey=st.get("standard_inchi_key"),
            formula=pr.get("full_molformula"),
            molecular_weight=_float(pr.get("full_mwt")), monoisotopic_mass=_float(pr.get("mw_monoisotopic")),
            xlogp=_float(pr.get("alogp")), tpsa=_float(pr.get("psa")),
            hbd=_int(pr.get("hbd")), hba=_int(pr.get("hba")), rotatable_bonds=_int(pr.get("rtb")),
            extra={k: v for k, v in extra.items() if v is not None},
        )

    def by_inchikey(self, inchikey: str) -> Optional[CompoundRecord]:
        m = self._get(f"{_BASE}/molecule/{inchikey}.json")
        return self._to_record(m, inchikey) if m and m.get("molecule_chembl_id") else None

    def by_id(self, identifier: str) -> Optional[CompoundRecord]:
        m = self._get(f"{_BASE}/molecule/{identifier.upper()}.json")
        return self._to_record(m, identifier) if m and m.get("molecule_chembl_id") else None

    def by_name(self, name: str) -> Optional[CompoundRecord]:
        # exact pref_name first, then the full-text search
        data = self._get(f"{_BASE}/molecule.json", {"pref_name__iexact": name, "limit": 1})
        mols = (data or {}).get("molecules") or []
        if not mols:
            data = self._get(f"{_BASE}/molecule/search.json", {"q": name, "limit": 1})
            mols = (data or {}).get("molecules") or []
        return self._to_record(mols[0], name) if mols else None

    def by_smiles(self, smiles: str) -> Optional[CompoundRecord]:
        data = self._get(f"{_BASE}/molecule.json", {"molecule_structures__canonical_smiles__flexmatch": smiles, "limit": 1})
        mols = (data or {}).get("molecules") or []
        if mols:
            return self._to_record(mols[0], smiles)
        return super().by_smiles(smiles)
