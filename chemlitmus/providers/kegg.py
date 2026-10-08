"""KEGG COMPOUND provider (KEGG REST, flat-file format)."""

from __future__ import annotations

import re
from typing import Dict, List, Optional

from chemlitmus.providers.base import CompoundRecord, Provider, _float

_BASE = "https://rest.kegg.jp"
_ID_RE = re.compile(r"^(?:cpd:|dr:)?[CD]\d{5}$", re.I)


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


_QUALIFIER = re.compile(r"\s*\((?:JP\d+|USP|INN|JAN|BAN|TN|USAN|DCF|NF)(?:/(?:JP\d+|USP|INN|JAN|BAN|TN|USAN|DCF|NF))*\)\s*$")


def _strip_qualifier(n: str) -> str:
    return _QUALIFIER.sub("", n.strip())


def _structure_from_molblock(block: str):
    """SMILES, InChI and InChIKey from a KEGG MOL block (KEGG publishes no InChIKey itself)."""
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromMolBlock(block, sanitize=True)
        if mol is None:
            mol = Chem.MolFromMolBlock(block, sanitize=False)
            if mol is None:
                return None, None, None
            mol.UpdatePropertyCache(strict=False)
        smiles = Chem.MolToSmiles(mol)
        inchi = Chem.MolToInchi(mol) or None
        return smiles, inchi, (Chem.InchiToInchiKey(inchi) if inchi else None)
    except Exception:
        return None, None, None


class KEGGProvider(Provider):
    key = "kegg"
    name = "KEGG"
    homepage = "https://www.kegg.jp"
    min_interval = 0.4                                     # KEGG asks for ≤3 requests/s

    def looks_like_id(self, query: str) -> bool:
        return bool(_ID_RE.match(query))

    def by_id(self, identifier: str) -> Optional[CompoundRecord]:
        cid = identifier.split(":")[-1].upper()
        if cid.startswith("D") and self._drug_same_as(cid):
            cid = self._drug_same_as(cid)  # a DRUG entry that is literally a COMPOUND entry
        text = self._get(f"{_BASE}/get/{cid}", json=False)
        if not text:
            return None
        f = _parse_flat(text)
        names = [_strip_qualifier(n.rstrip(";")) for n in f.get("NAME", [])]
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
        smiles = inchi = inchikey = None
        mol_block = self._get(f"{_BASE}/get/{cid}/mol", json=False)
        if mol_block and mol_block.strip():
            smiles, inchi, inchikey = _structure_from_molblock(mol_block)
        return CompoundRecord(
            source=self.key, source_id=cid, query=identifier,
            url=f"{self.homepage}/entry/{cid}",
            name=names[0] if names else None, synonyms=names[1:25],
            smiles=smiles, inchi=inchi, inchikey=inchikey,
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

    def _drug_same_as(self, did: str) -> Optional[str]:
        text = self._get(f"{_BASE}/get/{did}", json=False)
        if not text:
            return None
        for ln in _parse_flat(text).get("REMARK", []):
            m = re.search(r"Same as:\s*(C\d{5})", ln)
            if m:
                return m.group(1)
        return None

    def by_name(self, name: str) -> Optional[CompoundRecord]:
        """Exact (case-insensitive) name match in KEGG COMPOUND, then KEGG DRUG.

        ``find`` is a substring search — "risperidone" returns *9-hydroxyrisperidone* — so a hit
        whose name list does not contain the query verbatim is not accepted; ``resolve`` will then
        recover the entry through the consensus structure instead. KEGG DRUG names carry
        qualifiers such as "(JP19/USP/INN)" and "(TN)", which are ignored for matching.
        """
        lname = name.strip().lower()
        for db in ("compound", "drug"):
            text = self._get(f"{_BASE}/find/{db}/{name}", json=False)
            if not text or not text.strip():
                continue
            for line in text.splitlines():
                kid, _, names = line.partition("\t")
                cands = {_strip_qualifier(n).lower() for n in names.split(";")}
                if lname in cands:
                    return self.by_id(kid.split(":")[-1])
        return None

    def by_inchikey(self, inchikey: str) -> Optional[CompoundRecord]:
        # KEGG has no InChIKey endpoint; resolve through ChEBI's cross-references when available.
        try:
            from chemlitmus.providers.chebi import ChEBIProvider
            ch = ChEBIProvider().by_inchikey(inchikey)
        except Exception:
            ch = None
        kid = (ch.cross_refs.get("kegg") if ch else None)
        return self.by_id(kid) if kid else None
