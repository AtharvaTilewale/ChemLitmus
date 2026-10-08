"""PubChem PUG REST API client for ChemLitmus."""

import time
import threading
import urllib.parse
from typing import Optional, Dict, Any, Union
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from pydantic import BaseModel, Field

from chemlitmus.logging_config import logger
from chemlitmus.config import settings

# -------------------------------------------------------------------------
# Data Models
# -------------------------------------------------------------------------

class PubChemCompound(BaseModel):
    """Data model representing a fetched compound's metadata."""
    input_query: str = Field(description="The original query string used to find this compound")
    cid: Optional[int] = None
    title: Optional[str] = None
    iupac_name: Optional[str] = None
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    canonical_smiles: Optional[str] = None
    isomeric_smiles: Optional[str] = None
    inchikey: Optional[str] = None
    inchi: Optional[str] = None
    xlogp: Optional[float] = None
    hbond_donor_count: Optional[int] = None
    hbond_acceptor_count: Optional[int] = None

# -------------------------------------------------------------------------
# PubChem API Client
# -------------------------------------------------------------------------

class PubChemClient:
    """
    Client for interacting with the PubChem PUG REST API.
    Includes thread-safe rate limiting and connection retries.
    """
    
    # PubChem allows max 5 requests per second. 
    # We use a thread lock to ensure batch multithreading respects this limit.
    _rate_limit_lock = threading.Lock()
    _last_request_time = 0.0
    _MIN_REQUEST_INTERVAL = 0.22  # slightly over 200ms to be safe

    def __init__(self):
        self.base_url = settings.pubchem_base_url.rstrip("/")
        
        # Define properties we want to fetch for every compound
        self.properties = (
            "Title,MolecularFormula,MolecularWeight,CanonicalSMILES,"
            "IsomericSMILES,ConnectivitySMILES,SMILES,InChI,InChIKey,IUPACName,XLogP,"
            "HBondDonorCount,HBondAcceptorCount"
        )
        
        # Setup robust session with Retry logic for transient server errors (500, 502, 503, 504)
        self.session = requests.Session()
        retries = Retry(
            total=3,
            backoff_factor=1,  # 1s, 2s, 4s delays between retries
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET"]
        )
        adapter = HTTPAdapter(max_retries=retries)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

    def _enforce_rate_limit(self):
        """Thread-safe rate limiter to prevent API blocks."""
        with self._rate_limit_lock:
            current_time = time.time()
            elapsed = current_time - self.__class__._last_request_time
            if elapsed < self._MIN_REQUEST_INTERVAL:
                time.sleep(self._MIN_REQUEST_INTERVAL - elapsed)
            self.__class__._last_request_time = time.time()

    def _make_request(self, url: str) -> Optional[Dict[str, Any]]:
        """Execute the HTTP GET request with rate limiting and error handling."""
        self._enforce_rate_limit()
        
        try:
            logger.debug(f"Requesting PubChem URL: {url}")
            response = self.session.get(url, timeout=10)
            
            if response.status_code == 200:
                return response.json()
            elif response.status_code == 404:
                logger.debug(f"PubChem returned 404 Not Found for URL: {url}")
                return None
            else:
                logger.warning(f"PubChem API error {response.status_code}: {response.text}")
                return None
                
        except requests.exceptions.RequestException as e:
            logger.error(f"Network error querying PubChem: {e}")
            return None

    def _parse_response(self, data: Dict[str, Any], input_query: str) -> Optional[PubChemCompound]:
        """Safely parse the PubChem JSON response into a Pydantic model."""
        try:
            properties = data.get("PropertyTable", {}).get("Properties", [])
            if not properties:
                return None
                
            # Usually we only care about the first (best) match returned
            prop = properties[0]
            
            canonical_smiles = (
                prop.get("CanonicalSMILES")
                or prop.get("SMILES")
                or prop.get("ConnectivitySMILES")
            )
            isomeric_smiles = (
                prop.get("IsomericSMILES")
                or prop.get("SMILES")
                or prop.get("CanonicalSMILES")
                or prop.get("ConnectivitySMILES")
            )

            return PubChemCompound(
                input_query=str(input_query),
                cid=prop.get("CID"),
                title=prop.get("Title"),
                iupac_name=prop.get("IUPACName"),
                molecular_formula=prop.get("MolecularFormula"),
                molecular_weight=prop.get("MolecularWeight"),
                canonical_smiles=canonical_smiles,
                isomeric_smiles=isomeric_smiles,
                inchikey=prop.get("InChIKey"),
                inchi=prop.get("InChI"),
                xlogp=prop.get("XLogP"),
                hbond_donor_count=prop.get("HBondDonorCount"),
                hbond_acceptor_count=prop.get("HBondAcceptorCount")
            )
        except Exception as e:
            logger.error(f"Failed to parse PubChem response for {input_query}: {e}")
            return None

    # -------------------------------------------------------------------------
    # Specific Lookup Methods
    # -------------------------------------------------------------------------

    def lookup_by_cid(self, cid: Union[int, str]) -> Optional[PubChemCompound]:
        """Lookup compound by PubChem CID."""
        url = f"{self.base_url}/compound/cid/{cid}/property/{self.properties}/JSON"
        data = self._make_request(url)
        if data:
            return self._parse_response(data, input_query=str(cid))
        return None

    def lookup_by_smiles(self, smiles: str) -> Optional[PubChemCompound]:
        """Lookup compound by SMILES string."""
        encoded_smiles = urllib.parse.quote(smiles)
        url = f"{self.base_url}/compound/smiles/{encoded_smiles}/property/{self.properties}/JSON"
        data = self._make_request(url)
        if data:
            return self._parse_response(data, input_query=smiles)
        return None

    def lookup_by_name(self, name: str) -> Optional[PubChemCompound]:
        """Lookup compound by common or IUPAC name."""
        encoded_name = urllib.parse.quote(name)
        url = f"{self.base_url}/compound/name/{encoded_name}/property/{self.properties}/JSON"
        data = self._make_request(url)
        if data:
            return self._parse_response(data, input_query=name)
        return None

    def lookup_by_inchikey(self, inchikey: str) -> Optional[PubChemCompound]:
        """Lookup compound by InChIKey."""
        url = f"{self.base_url}/compound/inchikey/{inchikey}/property/{self.properties}/JSON"
        data = self._make_request(url)
        if data:
            return self._parse_response(data, input_query=inchikey)
        return None

    # -------------------------------------------------------------------------
    # Main Router
    # -------------------------------------------------------------------------

    def lookup(self, query: str, search_type: str = "auto") -> Optional[PubChemCompound]:
        """
        Main router method. Automatically detects input type or uses specified type.
        
        Args:
            query: The string or number to search for.
            search_type: "auto", "cid", "smiles", "name", or "inchikey".
            
        Returns:
            PubChemCompound if found, else None.
        """
        query_str = str(query).strip()
        search_type = search_type.lower()
        
        logger.debug(f"Executing PubChem lookup for '{query_str}' (type: {search_type})")

        # 1. Explicit CID or all digits in auto mode
        if search_type == "cid" or (search_type == "auto" and query_str.isdigit()):
            return self.lookup_by_cid(query_str)
            
        # 2. Explicit searches
        if search_type == "name":
            return self.lookup_by_name(query_str)
        if search_type == "smiles":
            return self.lookup_by_smiles(query_str)
        if search_type == "inchikey":
            return self.lookup_by_inchikey(query_str)

        # 3. Smart Auto-Detection Mode
        if search_type == "auto":
            # If it's 27 characters and has two hyphens, it's likely an InChIKey
            if len(query_str) == 27 and query_str[14] == "-" and query_str[25] == "-":
                return self.lookup_by_inchikey(query_str)
                
            # Attempt to validate as SMILES via RDKit
            # We import here to avoid circular imports if pubchem.py is loaded early
            try:
                from chemlitmus.core.smiles import validate_smiles
                validation = validate_smiles(query_str)
                if validation.is_valid:
                    return self.lookup_by_smiles(validation.canonical_smiles)
            except ImportError:
                logger.warning("Could not import validate_smiles. Skipping SMILES validation.")
                pass
                
            # If it's not a digit, not an InChIKey, and not a valid SMILES...
            # The safest fallback is treating it as a chemical Name!
            return self.lookup_by_name(query_str)
            
        logger.error(f"Unknown search_type: {search_type}")
        return None
    
    def download_structure(self, cid: int, format: str = "sdf", dimension: str = "3d", output_dir: str = "structures", force: bool = False) -> str:
        """
        Download physical structure file from PubChem.
        Includes smart resume logic to skip existing files.
        """
        from pathlib import Path
        
        format = format.lower()
        dimension = dimension.lower()
        
        # Build the save path
        out_path = Path(output_dir) / f"{cid}_{dimension}.{format}"
        
        # Smart resume logic
        if out_path.exists() and not force:
            return "Skipped (File already exists)"
            
        # Construct the exact URL based on requested format and dimension
        if format == "png":
            # Images don't use the record_type parameter in the same way
            url = f"{self.base_url}/compound/cid/{cid}/PNG?image_size=large"
        else:
            url = f"{self.base_url}/compound/cid/{cid}/{format.upper()}"
            # PubChem requires the record_type query parameter for 3D structures
            url += "?record_type=3d" if dimension == "3d" else "?record_type=2d"
            
        self._enforce_rate_limit()
        
        try:
            # We use self.session to inherit the retry logic and headers
            response = self.session.get(url, timeout=15)
            
            if response.status_code == 200:
                out_path.parent.mkdir(parents=True, exist_ok=True)
                # Write in binary mode to support images (PNG) and text (SDF) safely
                with open(out_path, "wb") as f:
                    f.write(response.content)
                return "Downloaded"
                
            elif response.status_code == 404:
                return "Not Found (3D structure might not be computed for this CID yet)"
            elif response.status_code == 400:
                return f"Bad Request (Format '{format}' might not be supported)"
            else:
                return f"HTTP Error {response.status_code}"
                
        except requests.exceptions.RequestException as e:
            logger.error(f"Network error downloading CID {cid}: {e}")
            return f"Network Error: {e}"