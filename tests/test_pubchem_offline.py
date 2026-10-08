"""Offline tests for the PubChem client, result export and UniChem (HTTP mocked)."""

import json
from types import SimpleNamespace

import pytest
import requests

from chemlitmus.core.pubchem import PubChemClient, PubChemCompound
from chemlitmus.providers import unichem as uc
from chemlitmus.providers.base import ProviderError
from chemlitmus.utils.export import export_results

PROPS = {"PropertyTable": {"Properties": [{"CID": 241, "Title": "Benzene", "MolecularFormula": "C6H6", "MolecularWeight": "78.11",
         "SMILES": "C1=CC=CC=C1", "InChIKey": "UHOVQNZJYSORNB-UHFFFAOYSA-N", "InChI": "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H",
         "IUPACName": "benzene", "XLogP": 2.1, "HBondDonorCount": 0, "HBondAcceptorCount": 0}]}}


class FakeResponse:
    def __init__(self, status, payload=None, text=""):
        self.status_code, self._payload, self.text = status, payload, text
    def json(self):
        return self._payload


@pytest.fixture
def fake_session(monkeypatch):
    calls = []
    def fake_get(url, timeout=10, **kw):
        calls.append(url)
        if "/cid/404/" in url or "name/nothing/" in url:
            return FakeResponse(404, text="not found")
        if "name/busy/" in url:
            return FakeResponse(503, text="PUGREST.ServerBusy")
        if "name/boom/" in url:
            raise requests.exceptions.ConnectionError("down")
        return FakeResponse(200, PROPS)
    monkeypatch.setattr(PubChemClient, "_enforce_rate_limit", lambda self: None)
    client = PubChemClient()
    monkeypatch.setattr(client, "session", SimpleNamespace(get=fake_get))
    return client, calls


def test_pubchem_parse_and_routing(fake_session):
    client, calls = fake_session
    c = client.lookup("241", "auto")                       # digits -> CID
    assert isinstance(c, PubChemCompound) and c.cid == 241 and c.canonical_smiles == "C1=CC=CC=C1"
    assert c.isomeric_smiles == "C1=CC=CC=C1" and c.iupac_name == "benzene" and "/cid/241/" in calls[-1]
    assert client.lookup("c1ccccc1").cid == 241 and "/smiles/" in calls[-1]                    # valid SMILES
    assert client.lookup("UHOVQNZJYSORNB-UHFFFAOYSA-N").cid == 241 and "/inchikey/" in calls[-1]
    assert client.lookup("benzene").cid == 241 and "/name/" in calls[-1]                       # fallback: name
    assert client.lookup("benzene", "name").cid == 241
    assert client.lookup("c1ccccc1", "smiles").cid == 241
    assert client.lookup("UHOVQNZJYSORNB-UHFFFAOYSA-N", "inchikey").cid == 241
    assert client.lookup("x", "bogus") is None


def test_pubchem_error_paths(fake_session):
    client, _ = fake_session
    assert client.lookup("404", "cid") is None
    assert client.lookup("nothing", "name") is None
    assert client.lookup("busy", "name") is None
    assert client.lookup("boom", "name") is None
    assert client._parse_response({"PropertyTable": {"Properties": []}}, "q") is None
    assert client._parse_response({"garbage": 1}, "q") is None


def test_export_results(tmp_path):
    recs = [PubChemCompound(input_query="benzene", cid=241, title="Benzene", molecular_formula="C6H6"),
            PubChemCompound(input_query="ethanol", cid=702, title="Ethanol", molecular_formula="C2H6O")]
    export_results(recs, tmp_path / "r.csv", "csv")
    assert (tmp_path / "r.csv").read_text().count("\n") == 3
    export_results(recs, tmp_path / "r.json", "JSON")
    assert json.loads((tmp_path / "r.json").read_text())[1]["cid"] == 702
    export_results(recs, tmp_path / "r.xlsx", "xlsx")
    assert (tmp_path / "r.xlsx").stat().st_size > 0
    export_results([], tmp_path / "empty.csv", "csv")
    assert not (tmp_path / "empty.csv").exists()
    with pytest.raises(ValueError):
        export_results(recs, tmp_path / "r.txt", "txt")


def test_unichem_xrefs_parsing(monkeypatch):
    payload = [{"src_id": "1", "src_compound_id": "CHEMBL25"}, {"src_id": "22", "src_compound_id": "2244"},
               {"src_id": "22", "src_compound_id": "9999"}, {"src_id": "999", "src_compound_id": "Z"}]
    monkeypatch.setattr(uc, "_last", 0.0)
    monkeypatch.setattr(uc.requests, "get", lambda url, timeout, **kw: FakeResponse(200, payload))
    x = uc.unichem_xrefs("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
    assert x["chembl"] == "CHEMBL25" and x["pubchem"] == "2244" and x["unichem_src_999"] == "Z"
    monkeypatch.setattr(uc.requests, "get", lambda url, timeout, **kw: FakeResponse(404, {"error": "not found"}))
    assert uc.unichem_xrefs("AAAAAAAAAAAAAA-BBBBBBBBBB-N") == {}


def test_unichem_errors(monkeypatch):
    monkeypatch.setattr(uc, "_last", 0.0)
    monkeypatch.setattr(uc.requests, "get", lambda url, timeout, **kw: FakeResponse(500, None, "oops"))
    with pytest.raises(ProviderError):
        uc.unichem_xrefs("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
    def boom(url, timeout, **kw):
        raise requests.exceptions.ConnectionError("down")
    monkeypatch.setattr(uc.requests, "get", boom)
    with pytest.raises(ProviderError):
        uc.unichem_xrefs("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
