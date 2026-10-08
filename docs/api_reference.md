# API Reference

This page details the core public functions and data models exposed by the `chemlitmus` package. All core functions can be imported directly from the top-level package.

```python
from chemlitmus import lookup, lookup_by_name, lookup_file, validate_smiles, download_structure, generate_structure
```

`lookup()`

The primary router for fetching chemical data. It auto-detects the format of the query and routes it to the appropriate PubChem endpoint.

- ### Parameters:
    - `query` (str | int): The identifier to search for (SMILES, CID, InChIKey, or Name).
    - `search_type` (str, optional): Override auto-detection. Valid options are "`auto`", "`cid`", "`smiles`", "`name`", or "`inchikey`". Defaults to "`auto`".
    - `use_cache` (bool, optional): Whether to check the local SQLite database before making a network request. Defaults to `True`.

- ### Returns:
    - `PubChemCompound` if found, otherwise `None`.

`lookup_by_name()`

```python
def lookup_by_name(name: str, use_cache: bool = True) -> Optional[PubChemCompound]
```

Explicitly queries PubChem by a common or IUPAC chemical name. This bypasses all SMILES validation checks.

- ### Parameters:
    - `name` (str): The chemical name (e.g., "`Aspirin`", "`Benzene`").
    - `use_cache` (bool, optional): Defaults to `True`.

- ### Returns:
    - `PubChemCompound` if found, otherwise `None`.

`lookup_file()`

```python
def lookup_file(input_file: Union[str, Path], output_file: Optional[Union[str, Path]] = None, output_format: str = "csv", remove_duplicates: bool = True) -> List[PubChemCompound]
```

Processes a batch file of compounds using multithreading.

- ### Parameters:
    - `input_file` (str | Path): Path to the input file (`.csv`, `.tsv`, `.xlsx`, `.smi`, `.sdf`).
    - `output_file` (str | Path, optional): Path to save the results. If `None`, results are kept in memory.
    - `output_format` (str, optional): Format to export ("`csv`", "`xlsx`", "`json`").
    - `remove_duplicates` (bool, optional): Automatically dedupes the input list to save API calls. Defaults to `True`.

- ### Returns:
    - A list of `PubChemCompound` objects.
  
`download_structure()`

```python
def download_structure(cid: int, format: str = "sdf", dimension: str = "3d", output_dir: str = "structures", force: bool = False) -> str
```

Downloads physical structure files directly from PubChem. Features smart-resume to skip existing files.

- ### Parameters:
    - `cid` (int): The PubChem CID.
    - `format` (str, optional): "`sdf`", "`mol`", "`pdb`", or "`png`". Defaults to "`sdf`".
    - `dimension` (str, optional): "`2d`" or "`3d`". Defaults to "`3d`".
    - `output_dir` (str, optional): The folder to save structures in. Defaults to "`structures`".
    - `force` (bool, optional): Overwrite existing files. Defaults to `False`.

- ### Returns:
    - A status string ("`Downloaded`", "`Skipped`", or error message).

`generate_structure()`

```python
def generate_structure(smiles: str, output_path: Union[str, Path], format: str = "sdf", dimension: str = "3d", force: bool = False, title: Optional[str] = None) -> str
```

Generates 2D or 3D molecular conformations offline from a SMILES string using RDKit and saves them to file.

- ### Parameters:
    - `smiles` (str): Input SMILES string.
    - `output_path` (str | Path): Output file destination path.
    - `format` (str, optional): Supported formats are "`sdf`", "`mol`", and "`pdb`". Defaults to "`sdf`".
    - `dimension` (str, optional): "`2d`" or "`3d`". Defaults to "`3d`".
    - `force` (bool, optional): Overwrite existing file if `True`. Defaults to `False`.
    - `title` (str, optional): Compound title or identifier to embed in the structure.

- ### Returns:
    - A status string ("`Generated`", "`Skipped (File already exists)`", or error message).

`validate_smiles()`

```python
def validate_smiles(smiles_str: str) -> SMILESValidationResult
```

Performs high-speed, offline SMILES validation and physicochemical descriptor calculation using the local RDKit engine. Does not connect to the internet.

- ### Parameters:
    - `smiles_str` (str): The SMILES string to validate.

- ### Returns:
    - A `SMILESValidationResult` object.

## Data Models

ChemLitmus uses [Pydantic](https://pydantic-docs.helpmanual.io/) models to strictly type and validate returned data.

`PubChemCompound`

The core object returned by all `lookup` methods.

| Attribute | Type | Description |
|-----------|------|-------------|
| `input_query` | `str` | The original string used to search for this compound. |
| `cid` | `int` | The official PubChem Compound ID. |
| `title` | `str` | The common title of the compound. |
| `iupac_name` | `str` | The IUPAC standardized name. |
| `molecular_formula` | `str` | The chemical formula (e.g., `C6H6`). |
| `molecular_weight` | `float` | Exact molecular weight in g/mol. |
| `canonical_smiles` | `str` | PubChem's standardized canonical SMILES. |
| `isomeric_smiles` | `str` | PubChem's standardized isomeric SMILES. |
| `inchikey` | `str` | The 27-character InChIKey. |
| `inchi` | `str` | The full InChI string. |
| `xlogp` | `float` | Calculated lipophilicity (XLogP3). |
| `hbond_donor_count` | `int` | Number of hydrogen bond donors. |
| `hbond_acceptor_count` | `int` | Number of hydrogen bond acceptors. |

`SMILESValidationResult`

The object returned by the offline `validate_smiles()` function.

| Attribute | Type | Description |
|-----------|------|-------------|
| `is_valid` | `bool` | `True` if RDKit successfully parsed the SMILES string. |
| `canonical_smiles` | `str` | RDKit-standardized canonical SMILES. |
| `logp` | `float` | RDKit-calculated MolLogP value. |
| `tpsa` | `float` | Topological Polar Surface Area (TPSA). |
| `heavy_atom_count` | `int` | Total number of heavy (non-hydrogen) atoms. |
| `error_message` | `str` | Error message returned if `is_valid` is `False`. |

---

## Cheminformatics Functions

All three functions below are **fully offline** and powered by RDKit. No network access is required.

```python
from chemlitmus import compute_fingerprint, apply_filters, compute_similarity
```

---

`compute_fingerprint()`

```python
def compute_fingerprint(smiles: str, fp_type: str = "ecfp4", n_bits: int = 2048) -> Union[FingerprintResult, List[FingerprintResult]]
```

Computes molecular fingerprints from a SMILES string. Uses the modern `rdFingerprintGenerator` API internally.

- ### Parameters:
    - `smiles` (str): Input SMILES string.
    - `fp_type` (str): Fingerprint algorithm. Options: `"ecfp4"`, `"ecfp6"`, `"fcfp4"`, `"maccs"`, `"rdkit"`, `"atompair"`, `"torsion"`, or `"all"`.
    - `n_bits` (int): Number of bits for hashed fingerprints. Ignored for MACCS (fixed at 167). Default `2048`.

- ### Returns:
    - A single `FingerprintResult`, or a `List[FingerprintResult]` when `fp_type="all"`.

---

`apply_filters()`

```python
def apply_filters(smiles: str, rules: Optional[List[str]] = None) -> FilterResult
```

Evaluates drug-likeness and ADMET filters on a SMILES string.

- ### Parameters:
    - `smiles` (str): Input SMILES string.
    - `rules` (list, optional): Rule names to apply. Valid values: `"lipinski"`, `"veber"`, `"ghose"`, `"egan"`, `"ro3"`, `"pains"`, `"qed"`. Pass `None` or `["all"]` to apply every rule.

- ### Returns:
    - A `FilterResult` with per-rule `RuleResult` objects and computed property values.

---

`compute_similarity()`

```python
def compute_similarity(query_smiles: str, library: List[str], fp_type: str = "ecfp4", n_bits: int = 2048, threshold: float = 0.0, top_n: Optional[int] = None) -> List[SimilarityResult]
```

Computes Tanimoto similarity between a query SMILES and a list of library SMILES.

- ### Parameters:
    - `query_smiles` (str): Query compound SMILES.
    - `library` (List[str]): List of library SMILES to compare against.
    - `fp_type` (str): Fingerprint type for comparison. Same options as `compute_fingerprint`.
    - `n_bits` (int): Fingerprint size in bits.
    - `threshold` (float): Minimum Tanimoto score to include (0.0–1.0). Default `0.0`.
    - `top_n` (int, optional): Return only the top N results. `None` returns all above threshold.

- ### Returns:
    - `List[SimilarityResult]` sorted by similarity descending.

---

## Data Models

`FingerprintResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `smiles` | `str` | Input SMILES string |
| `fingerprint_type` | `str` | Algorithm used (e.g., `"ecfp4"`) |
| `n_bits` | `int` | Total number of bits |
| `n_on_bits` | `int` | Number of set bits |
| `density` | `float` | Fraction of set bits (`n_on_bits / n_bits`) |
| `bit_string` | `str` | Full binary bit-string (0s and 1s) |
| `hex_string` | `str` | Hex-encoded fingerprint |

`FilterResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `smiles` | `str` | Input SMILES string |
| `molecular_weight` | `float` | Average molecular weight (g/mol) |
| `logp` | `float` | RDKit MolLogP |
| `hbd` | `int` | H-bond donor count |
| `hba` | `int` | H-bond acceptor count |
| `tpsa` | `float` | Topological Polar Surface Area (A^2) |
| `rotatable_bonds` | `int` | Rotatable bond count |
| `heavy_atom_count` | `int` | Heavy atom count |
| `molar_refractivity` | `float` | Molar refractivity (MR) |
| `qed_score` | `float` | QED drug-likeness score (0–1) |
| `lipinski` | `RuleResult` | Lipinski Ro5 pass/fail |
| `veber` | `RuleResult` | Veber oral bioavailability pass/fail |
| `ghose` | `RuleResult` | Ghose filter pass/fail |
| `egan` | `RuleResult` | Egan filter pass/fail |
| `ro3` | `RuleResult` | Rule of Three (lead-likeness) pass/fail |
| `pains` | `RuleResult` | PAINS alert detection result |
| `passes_all` | `bool` | `True` if all requested rules pass |
| `error` | `str` | Error message if SMILES is invalid |

`RuleResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `passed` | `bool` | `True` if compound satisfies this rule |
| `details` | `str` | Human-readable explanation |

`SimilarityResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `rank` | `int` | Rank (1 = most similar) |
| `query` | `str` | Query SMILES |
| `hit` | `str` | Library SMILES of the hit |
| `similarity` | `float` | Tanimoto similarity score (0–1) |
| `fingerprint_type` | `str` | Fingerprint algorithm used |
