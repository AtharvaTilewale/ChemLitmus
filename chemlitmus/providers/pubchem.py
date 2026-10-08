"""PubChem provider (PUG REST)."""

from __future__ import annotations

from typing import Any, Dict, Optional
from urllib.parse import quote

from chemlitmus.config import settings
from chemlitmus.providers.base import CompoundRecord, Provider, _float, _int

_PROPS = ("Title,MolecularFormula,MolecularWeight,MonoisotopicMass,Charge,CanonicalSMILES,IsomericSMILES,"
          "SMILES,InChI,InChIKey,IUPACName,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount")


class PubChemProvider(Provider):
    key = "pubchem"
    name = "PubChem"
    homepage = "https://pubchem.ncbi.nlm.nih.gov"
    min_interval = 0.22                                    # PubChem: ≤5 requests/s

    def __init__(self) -> None:
        self.retries = int(settings.pubchem_retries)
        super().__init__()
        self.base = settings.pubchem_base_url.rstrip("/")
        self.timeout = float(settings.pubchem_timeout)
        self.retries = int(settings.pubchem_retries)

    def looks_like_id(self, query: str) -> bool:
        return query.isdigit()

    def _fetch(self, path: str, query: str) -> Optional[CompoundRecord]:
        data = self._get(f"{self.base}/compound/{path}/property/{_PROPS}/JSON")
        if not data:
            return None
        props = (data.get("PropertyTable") or {}).get("Properties") or []
        if not props:
            return None
        return self._to_record(props[0], query)

    def _to_record(self, p: Dict[str, Any], query: str) -> CompoundRecord:
        cid = str(p.get("CID", ""))
        smiles = p.get("IsomericSMILES") or p.get("SMILES") or p.get("CanonicalSMILES")
        return CompoundRecord(
            source=self.key, source_id=cid, query=query,
            url=f"{self.homepage}/compound/{cid}" if cid else None,
            name=p.get("Title"), synonyms=[n for n in [p.get("IUPACName")] if n],
            smiles=smiles, inchi=p.get("InChI"), inchikey=p.get("InChIKey"),
            formula=p.get("MolecularFormula"),
            molecular_weight=_float(p.get("MolecularWeight")), monoisotopic_mass=_float(p.get("MonoisotopicMass")),
            charge=_int(p.get("Charge")),
            xlogp=_float(p.get("XLogP")), tpsa=_float(p.get("TPSA")),
            hbd=_int(p.get("HBondDonorCount")), hba=_int(p.get("HBondAcceptorCount")),
            rotatable_bonds=_int(p.get("RotatableBondCount")),
            extra={"iupac_name": p.get("IUPACName")},
        )

    def by_inchikey(self, inchikey: str) -> Optional[CompoundRecord]:
        return self._fetch(f"inchikey/{inchikey}", inchikey)

    def by_name(self, name: str) -> Optional[CompoundRecord]:
        return self._fetch(f"name/{quote(name, safe='')}", name)

    def by_id(self, identifier: str) -> Optional[CompoundRecord]:
        return self._fetch(f"cid/{identifier}", identifier)

    def by_smiles(self, smiles: str) -> Optional[CompoundRecord]:
        return self._fetch(f"smiles/{quote(smiles, safe='')}", smiles)
