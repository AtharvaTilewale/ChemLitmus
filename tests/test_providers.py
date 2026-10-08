"""Tests for the multi-database provider layer.

Offline tests mock HTTP at the Provider._get boundary; integration tests hit the live services.
"""

import json

import pytest
from typer.testing import CliRunner

from chemlitmus import PROVIDERS, DEFAULT_SOURCES, CompoundRecord, ResolveResult, get_provider, resolve
from chemlitmus.cli import app
from chemlitmus.providers.base import QUERY_TYPES, Provider, ProviderError, _looks_like_inchikey
from chemlitmus.providers.chebi import ChEBIProvider
from chemlitmus.providers.chembl import ChEMBLProvider
from chemlitmus.providers.kegg import KEGGProvider, _parse_flat
from chemlitmus.providers.pubchem import PubChemProvider
from chemlitmus.providers import unichem as uc

runner = CliRunner()
IK = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

# ----------------------------------------------------------------------------- fixtures: canned payloads

PUBCHEM = {"PropertyTable": {"Properties": [{"CID": 2244, "Title": "Aspirin", "MolecularFormula": "C9H8O4", "MolecularWeight": "180.16",
            "MonoisotopicMass": "180.04225873", "Charge": 0, "IsomericSMILES": "CC(=O)OC1=CC=CC=C1C(=O)O", "InChI": "InChI=1S/C9H8O4/...",
            "InChIKey": IK, "IUPACName": "2-acetyloxybenzoic acid", "XLogP": 1.2, "TPSA": 63.6, "HBondDonorCount": 1, "HBondAcceptorCount": 4, "RotatableBondCount": 3}]}}
CHEMBL = {"molecule_chembl_id": "CHEMBL25", "pref_name": "ASPIRIN", "max_phase": "4.0", "molecule_type": "Small molecule", "first_approval": 1950,
          "oral": True, "withdrawn_flag": False, "atc_classifications": ["B01AC06"],
          "molecule_structures": {"canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O", "standard_inchi": "InChI=1S/C9H8O4/...", "standard_inchi_key": IK},
          "molecule_properties": {"full_mwt": "180.16", "alogp": "1.31", "hba": 3, "hbd": 1, "psa": "63.60", "rtb": 2, "qed_weighted": "0.55", "num_ro5_violations": 0, "full_molformula": "C9H8O4"},
          "molecule_synonyms": [{"molecule_synonym": "Acetylsalicylic acid"}, {"molecule_synonym": "ASA"}]}
CHEBI = {"id": 15365, "chebi_accession": "CHEBI:15365", "name": "acetylsalicylic acid", "ascii_name": "acetylsalicylic acid", "stars": 3,
         "definition": "A member of the class of benzoic acids ...",
         "chemical_data": {"formula": "C9H8O4", "charge": 0, "mass": "180.159", "monoisotopic_mass": "180.04226"},
         "default_structure": {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "standard_inchi": "InChI=1S/C9H8O4/...", "standard_inchi_key": IK},
         "names": {"SYNONYM": [{"name": "aspirin"}, {"name": "ASA"}], "INN": [{"name": "aspirin"}], "BRAND NAME": [{"name": "Bayer"}]},
         "database_accessions": {"CAS": [{"source_name": "KEGG COMPOUND", "type": "CAS", "accession_number": "50-78-2"}],
                                 "MANUAL_X_REF": [{"source_name": "DrugBank", "type": "MANUAL_X_REF", "accession_number": "DB00945"},
                                                  {"source_name": "KEGG COMPOUND", "type": "MANUAL_X_REF", "accession_number": "C01405"},
                                                  {"source_name": "PDBeChem", "type": "MANUAL_X_REF", "accession_number": "AIN"}]}}
CHEBI_SEARCH_WRONG = {"results": [{"_score": 58.7, "_source": {"chebi_accession": "CHEBI:138615", "name": "aspirin-triggered resolvin D3", "stars": 3, "inchikey": "QBTJOLCUKWLTIC-AXRBAIMGSA-N"}},
                                  {"_score": 54.0, "_source": {"chebi_accession": "CHEBI:15365", "name": "acetylsalicylic acid", "stars": 3, "inchikey": IK}}]}
CHEBI_RESOLVIN = {"id": 138615, "chebi_accession": "CHEBI:138615", "name": "aspirin-triggered resolvin D3", "stars": 3,
                  "chemical_data": {"formula": "C22H32O5", "mass": "376.49"}, "default_structure": {"smiles": "CCCCC", "standard_inchi_key": "QBTJOLCUKWLTIC-AXRBAIMGSA-N"},
                  "names": {"SYNONYM": [{"name": "AT-RvD3"}]}, "database_accessions": {}}
KEGG_FLAT = """ENTRY       C01405                      Compound
NAME        Aspirin;
            Acetylsalicylic acid;
            2-Acetoxybenzenecarboxylic acid
FORMULA     C9H8O4
EXACT_MASS  180.0423
MOL_WEIGHT  180.1574
PATHWAY     map07048  Antimigraines
ENZYME      3.1.1.1
DBLINKS     CAS: 50-78-2
            PubChem: 4594
            ChEBI: 15365
            KNApSAcK: C00001492
///
"""
KEGG_MOL = """
  Marvin

 13 13  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2124    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4249    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6373    0.7000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4249   -1.4000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4249    2.8000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.6373    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8497    2.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6373    0.7000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  2  7  1  0  0  0  0
  7  8  2  0  0  0  0
  7  9  1  0  0  0  0
  3 10  1  0  0  0  0
 10 11  1  0  0  0  0
 11 12  1  0  0  0  0
 11 13  2  0  0  0  0
M  END
"""
KEGG_FIND = "cpd:C13400\tBufferin; Aspirin softam\ncpd:C01405\tAspirin; Acetylsalicylic acid; 2-Acetoxybenzenecarboxylic acid\n"
KEGG_DRUG_FIND = "dr:D00426\tRisperidone (JP19/USP/INN); Risperdal (TN)\n"
KEGG_DRUG_FLAT = "ENTRY       D00426                      Drug\nNAME        Risperidone (JP19/USP/INN);\nFORMULA     C23H27FN4O2\nDBLINKS     CAS: 106266-06-2\n            ChEBI: 8871\n///\n"
UNICHEM = [{"src_id": "1", "src_compound_id": "CHEMBL25"}, {"src_id": "2", "src_compound_id": "DB00945"}, {"src_id": "7", "src_compound_id": "CHEBI:15365"},
           {"src_id": "22", "src_compound_id": "2244"}, {"src_id": "84", "src_compound_id": "50-78-2"}, {"src_id": "999", "src_compound_id": "X1"}]


@pytest.fixture
def fake_http(monkeypatch):
    """Route every Provider._get call to a canned payload keyed on (provider, url fragment)."""
    calls = []
    def _get(self, url, params=None, *, json=True):
        calls.append((self.key, url, params))
        if self.key == "pubchem":
            if "/cid/404" in url or "name/nothing" in url:
                return None
            if "name/busy" in url:
                raise ProviderError("PubChem: HTTP 503 (PUGREST.ServerBusy)")
            return PUBCHEM
        if self.key == "chembl":
            if "/molecule/search.json" in url or ("molecule.json" in url and params and "pref_name__iexact" in params):
                return {"molecules": [CHEMBL]} if (params or {}).get("q", (params or {}).get("pref_name__iexact", "")).lower() in ("aspirin",) else {"molecules": []}
            if "flexmatch" in str(params):
                return {"molecules": [CHEMBL]}
            if url.endswith("/CHEMBL404.json"):
                return None
            return CHEMBL
        if self.key == "chebi":
            if "es_search" in url:
                term = (params or {}).get("term", "")
                return CHEBI_SEARCH_WRONG if term in ("aspirin", IK) else {"results": []}
            if "/compound/15365/" in url:
                return CHEBI
            if "/compound/138615/" in url:
                return CHEBI_RESOLVIN
            return None
        if self.key == "kegg":
            if "/find/compound/" in url:
                return KEGG_FIND if "aspirin" in url.lower() else ("cpd:C21516\tPaliperidone; 9-Hydroxyrisperidone\n" if "risperidone" in url.lower() else "")
            if "/find/drug/" in url:
                return KEGG_DRUG_FIND if "risperidone" in url.lower() else ""
            if url.endswith("/get/D00426"):
                return KEGG_DRUG_FLAT
            if url.endswith("/get/D00426/mol"):
                return ""
            if url.endswith("/get/C01405/mol"):
                return KEGG_MOL
            if "/get/C01405" in url:
                return KEGG_FLAT
            return None
        raise AssertionError(f"unexpected provider {self.key}")
    monkeypatch.setattr(Provider, "_get", _get)
    monkeypatch.setattr(uc, "unichem_xrefs", lambda ik, timeout=15.0: {uc.UNICHEM_SOURCES.get(r["src_id"], f"unichem_src_{r['src_id']}"): r["src_compound_id"] for r in UNICHEM})
    import chemlitmus.providers as prov
    monkeypatch.setattr(prov, "unichem_xrefs", uc.unichem_xrefs)
    return calls


# ----------------------------------------------------------------------------- unit: schema and helpers

def test_registry_and_constants():
    assert set(PROVIDERS) == {"pubchem", "chembl", "chebi", "kegg"}
    assert DEFAULT_SOURCES == ["pubchem", "chembl", "chebi"]
    assert QUERY_TYPES == ["auto", "name", "smiles", "inchikey", "id"]
    assert isinstance(get_provider("chembl"), ChEMBLProvider)
    assert get_provider("CHEMBL") is get_provider("chembl")          # cached, case-insensitive
    with pytest.raises(ValueError):
        get_provider("drugbank")


def test_looks_like_helpers():
    assert _looks_like_inchikey(IK) and not _looks_like_inchikey("CCO")
    assert PubChemProvider().looks_like_id("2244") and not PubChemProvider().looks_like_id("CHEMBL25")
    assert ChEMBLProvider().looks_like_id("chembl25") and not ChEMBLProvider().looks_like_id("2244")
    assert ChEBIProvider().looks_like_id("CHEBI:15365") and not ChEBIProvider().looks_like_id("15365")
    assert KEGGProvider().looks_like_id("C01405") and KEGGProvider().looks_like_id("cpd:C01405") and not KEGGProvider().looks_like_id("C1")


def test_kegg_flat_parser():
    f = _parse_flat(KEGG_FLAT)
    assert f["ENTRY"] == ["C01405                      Compound"]
    assert f["NAME"] == ["Aspirin;", "Acetylsalicylic acid;", "2-Acetoxybenzenecarboxylic acid"]
    assert f["DBLINKS"][2] == "ChEBI: 15365"


# ----------------------------------------------------------------------------- unit: each provider maps its payload

def test_pubchem_record(fake_http):
    r = PubChemProvider().lookup("2244")
    assert isinstance(r, CompoundRecord) and r.source == "pubchem" and r.source_id == "2244"
    assert r.name == "Aspirin" and r.inchikey == IK and r.molecular_weight == 180.16 and r.hba == 4
    assert r.url.endswith("/compound/2244") and r.extra["iupac_name"].startswith("2-acetyloxy")
    assert PubChemProvider().lookup("404", "id") is None


def test_chembl_record(fake_http):
    r = ChEMBLProvider().lookup("CHEMBL25")
    assert r.source_id == "CHEMBL25" and r.name == "ASPIRIN" and r.inchikey == IK
    assert r.extra["max_phase"] == "4.0" and r.extra["first_approval"] == 1950 and r.extra["qed_weighted"] == 0.55
    assert "ASA" in r.synonyms and r.formula == "C9H8O4" and r.tpsa == 63.6
    assert ChEMBLProvider().lookup("CHEMBL404", "id") is None
    assert ChEMBLProvider().lookup("CC(=O)Oc1ccccc1C(=O)O", "smiles").source_id == "CHEMBL25"


def test_chebi_record_and_xref_parsing(fake_http):
    r = ChEBIProvider().lookup("CHEBI:15365")
    assert r.source_id == "CHEBI:15365" and r.name == "acetylsalicylic acid" and r.charge == 0
    assert r.monoisotopic_mass == 180.04226 and r.extra["stars"] == 3
    # CAS entry is a registry number, not a KEGG id; MANUAL_X_REFs are real cross-references
    assert r.cross_refs == {"cas": "50-78-2", "drugbank": "DB00945", "kegg": "C01405", "pdbe": "AIN"}
    assert "aspirin" in r.synonyms and "Bayer" in r.synonyms


def test_chebi_name_search_prefers_synonym_match_over_score(fake_http):
    # top search hit is a derivative; the canonical entry lists 'aspirin' as a synonym
    r = ChEBIProvider().lookup("aspirin", "name")
    assert r.source_id == "CHEBI:15365"


def test_chebi_inchikey_requires_exact_key(fake_http):
    assert ChEBIProvider().lookup(IK).source_id == "CHEBI:15365"
    assert ChEBIProvider().lookup("AAAAAAAAAAAAAA-BBBBBBBBBB-N", "inchikey") is None


def test_kegg_record(fake_http):
    r = KEGGProvider().lookup("C01405")
    assert r.source_id == "C01405" and r.name == "Aspirin" and r.formula == "C9H8O4"
    assert r.monoisotopic_mass == 180.0423 and r.molecular_weight == 180.1574
    assert r.cross_refs == {"cas": "50-78-2", "pubchem_sid": "4594", "chebi": "CHEBI:15365", "knapsack": "C00001492"}
    assert r.extra["pathways"] == ["map07048"]
    assert r.inchikey == IK and r.smiles == "CC(=O)Oc1ccccc1C(=O)O"          # from the MOL block
    assert KEGGProvider().lookup("Aspirin", "name").source_id == "C01405"     # exact name, not the first substring hit
    assert KEGGProvider().lookup("nothing", "name") is None
    d = KEGGProvider().lookup("risperidone", "name")                          # COMPOUND has only a derivative; DRUG has it
    assert d.source_id == "D00426" and d.name == "Risperidone" and d.cross_refs["chebi"] == "CHEBI:8871"
    assert KEGGProvider().looks_like_id("D00426")


# ----------------------------------------------------------------------------- resolve()

def test_resolve_all_sources_agree(fake_http):
    res = resolve("aspirin", sources=["all"])
    assert isinstance(res, ResolveResult) and res.found
    assert {r.source for r in res.records} == {"pubchem", "chembl", "chebi", "kegg"}
    assert res.agreement == "agree" and res.consensus_inchikey == IK
    assert res.merged.name == "Aspirin" and res.merged.formula == "C9H8O4" and res.merged.inchikey == IK
    assert res.cross_refs["drugbank"] == "DB00945" and res.cross_refs["cas"] == "50-78-2" and res.cross_refs["kegg"] == "C01405"
    assert res.cross_refs["unichem_src_999"] == "X1"               # unknown UniChem sources are kept, not dropped
    assert res.by_source("chembl").extra["max_phase"] == "4.0"
    assert [o.source for o in res.outcomes] == ["pubchem", "chembl", "chebi", "kegg"]


def test_resolve_outlier_is_reresolved_by_structure(fake_http, monkeypatch):
    # make ChEBI's name search return only the wrong entry: resolve() must notice the InChIKey disagrees
    monkeypatch.setitem(CHEBI_SEARCH_WRONG, "results", CHEBI_SEARCH_WRONG["results"][:1])
    monkeypatch.setattr(ChEBIProvider, "by_inchikey", lambda self, ik: ChEBIProvider().by_id("CHEBI:15365") if ik == IK else None)
    res = resolve("aspirin", sources=["chembl", "chebi"], unichem=False)
    assert res.agreement == "agree" and res.by_source("chebi").source_id == "CHEBI:15365"


def test_resolve_not_found_and_errors(fake_http):
    res = resolve("nothing", sources=["pubchem", "chembl", "chebi", "kegg"], query_type="name", unichem=False)
    assert not res.found and res.merged is None and res.agreement == "unknown"
    assert all(o.status == "not found" for o in res.outcomes)
    res2 = resolve("busy", sources=["pubchem"], query_type="name", unichem=False)
    assert res2.outcomes[0].status == "error" and "503" in res2.outcomes[0].error


def test_resolve_validation():
    with pytest.raises(ValueError):
        resolve("x", sources=["drugbank"])
    with pytest.raises(ValueError):
        resolve("x", query_type="cas")


def test_resolve_json_roundtrip(fake_http):
    res = resolve("CHEMBL25", sources=["chembl"], unichem=False)
    data = json.loads(res.model_dump_json())
    assert ResolveResult.model_validate(data).records[0].source_id == "CHEMBL25"


def test_provider_cache_hits_by_query_id_and_inchikey(fake_http):
    p = ChEMBLProvider()
    assert p.lookup("aspirin", "name").source_id == "CHEMBL25"
    n_calls = len(fake_http)
    assert p.lookup("aspirin", "name").source_id == "CHEMBL25"          # same query: cached
    assert p.lookup("CHEMBL25").source_id == "CHEMBL25"                 # by id: cached via record
    assert p.lookup(IK).source_id == "CHEMBL25"                         # by inchikey: cached via record
    assert len(fake_http) == n_calls
    assert p.lookup("aspirin", "name", use_cache=False).source_id == "CHEMBL25"
    assert len(fake_http) == n_calls + 1
    from chemlitmus.core.database import DatabaseManager
    assert DatabaseManager().provider_cache_stats() == {"chembl": 3}


def test_resolve_second_run_is_offline(fake_http):
    resolve("aspirin", sources=["chembl", "chebi", "kegg"], unichem=False)
    n = len(fake_http)
    res = resolve("aspirin", sources=["chembl", "chebi", "kegg"], unichem=False)
    assert len(fake_http) == n and res.agreement == "agree"


# ----------------------------------------------------------------------------- disagreement classification

def test_resolve_explains_salt_form_disagreement(fake_http, monkeypatch):
    # ChEMBL answers with the sodium salt; ChEBI with the acid, whose synonyms include the query, so
    # ChEBI wins the tie and ChEMBL is re-queried by structure. Make that fail so the disagreement stands.
    rec_salt = ChEMBLProvider().lookup("CHEMBL25").model_copy(update={
        "source_id": "CHEMBL2260549", "name": "ASPIRIN SODIUM", "synonyms": [],
        "smiles": "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]", "inchikey": "SALTSALTSALTSA-UHFFFAOYSA-M"})
    monkeypatch.setattr(ChEMBLProvider, "by_name", lambda self, q: rec_salt)
    monkeypatch.setattr(ChEMBLProvider, "by_inchikey", lambda self, ik: None)
    res = resolve("aspirin", sources=["chembl", "chebi"], query_type="name", unichem=False, use_cache=False)
    assert res.agreement == "disagree" and len(res.records) == 2          # the uncorrected text hit is kept
    assert res.consensus_inchikey == IK
    assert res.agreement_level == "parent" and "salt" in res.disagreement
    assert res.pairwise[0].level == "parent"
    c = res.corrections[0]
    assert c.source == "chembl" and not c.corrected and "salt" in c.difference
    assert "disagrees" in res.outcomes[0].error


def test_resolve_records_successful_correction(fake_http, monkeypatch):
    monkeypatch.setitem(CHEBI_SEARCH_WRONG, "results", CHEBI_SEARCH_WRONG["results"][:1])
    monkeypatch.setattr(ChEBIProvider, "by_inchikey", lambda self, ik: ChEBIProvider().by_id("CHEBI:15365") if ik == IK else None)
    res = resolve("aspirin", sources=["chembl", "chebi"], unichem=False, use_cache=False)
    assert res.agreement == "agree" and res.agreement_level == "exact" and res.disagreement is None
    c = res.corrections[0]
    assert c.source == "chebi" and c.text_hit_id == "CHEBI:138615" and c.corrected and c.corrected_id == "CHEBI:15365"
    assert c.difference == "different compounds"


def test_agreement_level_none_vs_unavailable(fake_http, monkeypatch):
    res = resolve("CHEMBL25", sources=["chembl"], unichem=False, use_cache=False)
    assert res.agreement_level is None                                      # one structure: nothing to compare
    other = ChEMBLProvider().lookup("CHEMBL25").model_copy(update={"smiles": "CCCCCC", "inchikey": "WRONGWRONGWRON-UHFFFAOYSA-N"})
    monkeypatch.setattr(ChEBIProvider, "by_id", lambda self, q: other.model_copy(update={"source": "chebi", "source_id": "CHEBI:1"}))
    monkeypatch.setattr(ChEBIProvider, "by_inchikey", lambda self, ik: None)
    res = resolve("CHEMBL25", sources=["chembl", "chebi"], query_type="id", unichem=False, use_cache=False)
    assert res.agreement == "disagree" and res.agreement_level == "none" and res.disagreement == "different compounds"


def test_tie_break_prefers_source_that_knows_the_name(fake_http, monkeypatch):
    # ChEMBL's text search ranks a derivative first; ChEBI's record is literally named 'aspirin'.
    wrong = ChEMBLProvider().lookup("CHEMBL25").model_copy(update={"source_id": "CHEMBL999", "name": "ASPIRIN LYSINE",
                                                                   "synonyms": [], "smiles": "CCCCCC", "inchikey": "WRONGWRONGWRON-UHFFFAOYSA-N"})
    monkeypatch.setattr(ChEMBLProvider, "by_name", lambda self, q: wrong)
    monkeypatch.setattr(ChEMBLProvider, "by_inchikey", lambda self, ik: ChEMBLProvider().by_id("CHEMBL25") if ik == IK else None)
    res = resolve("aspirin", sources=["chembl", "chebi"], query_type="name", unichem=False, use_cache=False)
    assert res.consensus_inchikey == IK and res.agreement == "agree"
    assert res.corrections[0].source == "chembl" and res.corrections[0].corrected_id == "CHEMBL25"


def test_concordance_report(fake_http, monkeypatch):
    from chemlitmus.providers import concordance
    # "CHEMBL25" is unknown to ChEBI's text search but is recovered through the consensus InChIKey
    rep = concordance(["aspirin", "CHEMBL25", "nothing"], sources=["chembl", "chebi"], use_cache=False)
    assert rep.n_queries == 3 and rep.n_found_any == 2 and rep.n_found_all == 2
    assert rep.agreement_counts == {"agree": 2, "not found": 1}
    assert rep.level_counts == {"exact": 2}
    assert rep.found_by_source == {"chembl": 2, "chebi": 2}
    assert rep.rows[0].ids == {"chembl": "CHEMBL25", "chebi": "CHEBI:15365"}
    assert rep.rows[2].errors == {"chembl": "not found", "chebi": "not found"}


def test_cli_concordance(fake_http, tmp_path):
    f = tmp_path / "q.txt"; f.write_text("aspirin\nCHEMBL25\n")
    out, js = tmp_path / "c.csv", tmp_path / "c.json"
    r = runner.invoke(app, ["concordance", "--file", str(f), "-s", "chembl,chebi", "-o", str(out), "--json", str(js), "--no-cache"])
    assert r.exit_code == 0, r.output
    assert "Concordance across chembl, chebi" in r.output and "exact" in r.output
    rows = out.read_text().splitlines()
    assert rows[0].startswith("query,n_found,agreement,agreement_level") and len(rows) == 3
    assert json.loads(js.read_text())["n_queries"] == 2
    assert runner.invoke(app, ["concordance", "--file", str(f), "-s", "drugbank"]).exit_code == 1


# ----------------------------------------------------------------------------- CLI

def test_cli_resolve_single(fake_http, tmp_path):
    js = tmp_path / "r.json"
    r = runner.invoke(app, ["resolve", "aspirin", "--sources", "all", "--json", str(js)])
    assert r.exit_code == 0, r.output
    assert "sources agree" in r.output and "CHEMBL25" in r.output and "drugbank" in r.output
    assert json.loads(js.read_text())["agreement"] == "agree"


def test_cli_resolve_batch_and_errors(fake_http, tmp_path):
    f = tmp_path / "q.txt"; f.write_text("aspirin\nCHEMBL25\n")
    out = tmp_path / "r.csv"
    r = runner.invoke(app, ["resolve", "--file", str(f), "--sources", "chembl,chebi", "--no-unichem", "--output", str(out)])
    assert r.exit_code == 0, r.output
    rows = out.read_text().splitlines()
    assert rows[0].startswith("query,source,status") and len(rows) == 5         # header + 2 queries x 2 sources
    assert runner.invoke(app, ["resolve", "x", "--sources", "drugbank"]).exit_code == 1
    assert runner.invoke(app, ["resolve", "x", "--type", "cas"]).exit_code == 1
    assert runner.invoke(app, ["resolve"]).exit_code == 1
    assert runner.invoke(app, ["resolve", "nothing", "--type", "name", "--sources", "chembl"]).exit_code == 1


# ----------------------------------------------------------------------------- integration (network)

@pytest.mark.integration
def test_live_chembl_chebi_kegg_agree_on_aspirin():
    res = resolve("aspirin", sources=["chembl", "chebi", "kegg"], timeout=60)
    assert res.found and res.agreement == "agree" and res.consensus_inchikey == IK
    assert res.cross_refs.get("drugbank") == "DB00945"


@pytest.mark.integration
def test_live_unichem():
    x = uc.unichem_xrefs(IK)
    assert x.get("chembl") == "CHEMBL25" and x.get("chebi") == "CHEBI:15365"
