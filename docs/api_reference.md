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


## SMARTS Pattern Auditing

### `audit_smarts(patterns, library=None, library_source="user-supplied", checks=None, breadth_threshold=0.10, preparations=None, dead_sample=2500)`

Audit a set of SMARTS patterns against a reference molecule population.

| Parameter | Type | Description |
|-----------|------|-------------|
| `patterns` | `list[str]` or `list[tuple]` | SMARTS strings, or `(smarts, name, rule_set)` tuples from `load_patterns()` |
| `library` | `list[rdkit.Chem.Mol]` | Reference molecules; `None` loads the bundled set |
| `checks` | `list[str]` | Subset of `AUDIT_CHECKS` (`compile`, `breadth`, `dead`, `redundancy`, `sensitivity`); `None` runs all |
| `breadth_threshold` | `float` | Hit fraction above which a pattern is flagged over-broad |
| `preparations` | `list[str]` | Subset of `PREPARATIONS` for the sensitivity check; the first entry is the default preparation |
| `dead_sample` | `int` | Molecules sampled for atom-level triage of dead patterns |

Returns `SmartsAuditResult`.

### `explain_smarts(smarts, library=None, preparations=None, n_examples=5, sample=2500)`

Decompose one pattern: per-atom realisability, hits per preparation, example matches. Returns `SmartsExplanation`.

### `load_patterns(path)`

Read `(smarts, name, rule_set)` tuples from CSV/TSV/XLSX (column `smarts`, optional `name`/`description` and `rule_set_name`/`rule_set`) or from a text file with one SMARTS per line (optional trailing name; `#` comments).

### `load_reference_library(path=None, max_molecules=None)`

Returns `(molecules, source_label)`. `None` loads the bundled ChEMBL-derived set (9,272 molecules; see `chemlitmus/data/README.md`).

### `prepare_molecule(mol, preparation)`

Return a copy of `mol` as `implicit-h` (unchanged), `explicit-h` (`AddHs`) or `kekule` (Kekulé bonds, aromatic flags cleared).

### Data Models

`SmartsAuditResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_patterns`, `n_molecules` | `int` | Sizes |
| `library_source` | `str` | Where the reference molecules came from |
| `checks_run` | `list[str]` | Checks executed |
| `patterns` | `list[PatternAudit]` | One record per input pattern, in input order |
| `sensitivity` | `SensitivitySummary` | Catalogue-level reproducibility summary (`None` if not run) |
| `n_unparseable`, `n_needs_explicit_h`, `n_over_broad`, `n_dead`, `n_dead_never_matching_atom`, `n_duplicates`, `n_equivalent`, `n_subsumed`, `n_clean` | `int` | Aggregate counts (properties) |
| `to_rows()` | `list[dict]` | Flat rows for CSV export |

`PatternAudit`

| Attribute | Type | Description |
|-----------|------|-------------|
| `index`, `smarts`, `name`, `rule_set` | | Identity |
| `parses`, `parse_error` | `bool`, `str` | Compile status |
| `n_query_atoms`, `has_recursive_smarts` | `int`, `bool` | Structure |
| `requires_explicit_h` | `bool` | Contains a hydrogen query atom (positive, not negated or alternated) |
| `n_hits`, `hit_fraction`, `over_broad` | | Breadth under the default preparation |
| `never_fires`, `dead_verdict`, `never_matching_atoms` | | Dead-rule triage: `rare combination`, `never-matching atom`, or `fires only with <prep>` |
| `duplicate_of`, `equivalent_to`, `subsumed_by` | | Redundancy, as indices into `patterns` |
| `hits_by_preparation`, `preparation_sensitive` | `dict`, `bool` | Sensitivity |
| `flags`, `clean` | `list[str]`, `bool` | Summary (properties) |

`SensitivitySummary`

| Attribute | Type | Description |
|-----------|------|-------------|
| `preparations` | `list[str]` | Preparations compared |
| `compounds_flagged`, `total_hits`, `patterns_firing` | `dict[str,int]` | Per preparation |
| `verdict_flips` | `dict[str,int]` | Molecules whose pass/fail verdict differs from the default preparation |
| `n_sensitive_patterns` | `int` | Patterns whose hit count varies |

`SmartsExplanation`

| Attribute | Type | Description |
|-----------|------|-------------|
| `smarts`, `normalized_smarts`, `parses`, `parse_error` | | Identity and compile status |
| `n_query_atoms`, `n_query_bonds`, `requires_explicit_h`, `has_recursive_smarts` | | Structure |
| `atoms` | `list[AtomExplanation]` | Per atom: `query`, `n_matching_molecules`, `is_hydrogen` |
| `hits_by_preparation`, `n_molecules`, `example_matches` | | Whole-pattern behaviour |
| `never_matching_atoms`, `verdict` | | Diagnosis |

`FilterResult.preparation` and `SubstructureHit.preparation` record the molecule preparation used (`apply_filters(..., preparation=)`, `substructure_search(..., preparation=)`).


## Molecular Identity and Library Comparison

### `compute_identity(smiles)` → `IdentityKeys`

Identity keys at every level of `IDENTITY_LEVELS` (`exact`, `parent`, `tautomer`, `nostereo`, `skeleton`, `formula`). `parent` is the largest fragment after neutralisation; `tautomer`, `nostereo` and `skeleton` are RDKit `RegistrationHash` layers of the parent.

| Attribute | Type | Description |
|-----------|------|-------------|
| `input_smiles`, `is_valid`, `error` | | Record status |
| `exact`, `parent`, `tautomer`, `nostereo`, `skeleton`, `formula` | `str` | Keys |
| `n_fragments`, `had_charge`, `has_stereo` | | Properties of the input |
| `key(level)` | `str` | Key at a named level |

### `group_by_identity(smiles, level="parent", keys=None)` → `IdentityReport`

| Attribute | Type | Description |
|-----------|------|-------------|
| `level`, `n_records`, `n_valid` | | Inputs |
| `n_groups` | `int` | Distinct compounds at this level |
| `n_collapsed` | `int` | Records sharing a key with an earlier record |
| `n_groups_by_level` | `dict` | Distinct keys at every level |
| `groups` | `list[IdentityGroup]` | Multi-member groups (`key`, `size`, `indices`, `smiles`, `distinct_exact`, `differs_by`), largest first |
| `keys` | `list[IdentityKeys]` | Per-record keys |

### `strictest_shared_level(a, b)` / `describe_difference(a, b)`

For two `IdentityKeys`: the strictest level at which they agree (or `None`), and a plain-language description — `identical`, `salt, counter-ion or charge form`, `tautomer`, `stereochemistry`, `stereochemistry and tautomer`, `constitution (same formula only)`, `different compounds`.

### `diff_libraries(smiles_a, smiles_b, level="parent", keys_a=None, keys_b=None, include_unchanged=False)` → `LibraryDiff`

| Attribute | Type | Description |
|-----------|------|-------------|
| `n_a`, `n_b`, `n_valid_a`, `n_valid_b`, `n_keys_a`, `n_keys_b` | `int` | Sizes; `n_keys_*` are distinct compounds at the level |
| `n_added`, `n_removed`, `n_unchanged`, `n_changed` | `int` | Outcome counts |
| `changes_by_kind` | `dict[str,int]` | Changed compounds by reason |
| `multiplicity_changes` | `int` | Shared compounds whose record count differs |
| `jaccard` | `float` | Shared / all compounds (property) |
| `entries` | `list[DiffEntry]` | `status`, `key`, `smiles_a`, `smiles_b`, `count_a`, `count_b`, `change`; removed, then added, then changed |

## SMILES Diagnosis

### `diagnose_smiles(smiles, try_repair=True)` → `SmilesDiagnosis`

| Attribute | Type | Description |
|-----------|------|-------------|
| `input_smiles`, `is_valid`, `canonical_smiles` | | Status; a record with internal whitespace is reported invalid even though RDKit parses its first token |
| `problems` | `list[SmilesProblem]` | Ordered by `DIAGNOSTIC_CATEGORIES` then position |
| `repaired_smiles`, `repaired_is_valid`, `repairs_applied` | | Mechanical repair outcome, when attempted |
| `primary_category` | `str` | Category of the first problem (property) |
| `caret_line()` | `str` | `^` markers aligned under the input |

`SmilesProblem`: `category`, `message`, `position` (0-based), `length`, `atom_index`, `atom_indices`, `suggestion`.

`DIAGNOSTIC_CATEGORIES = ["characters", "brackets", "parentheses", "rings", "syntax", "valence", "aromaticity"]`.
