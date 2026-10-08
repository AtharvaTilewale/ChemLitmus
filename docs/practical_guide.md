# ChemLitmus Practical Guide

Welcome to the comprehensive tutorial for **ChemLitmus**. This guide is designed to take you from basic lookups to advanced, multithreaded batch processing using both the Command Line Interface (CLI) and the Python API.

---

# Table of Contents

- [ChemLitmus Practical Guide](#chemlitmus-practical-guide)
- [Table of Contents](#table-of-contents)
- [1. System Management](#1-system-management)
    - [CLI Commands](#cli-commands)
- [2. Single Compound Lookups](#2-single-compound-lookups)
  - [CLI Examples](#cli-examples)
    - [Basic Lookups:](#basic-lookups)
    - [Advanced CLI Flags:](#advanced-cli-flags)
    - [Python API Examples](#python-api-examples)
- [3. Batch Processing \& File I/O](#3-batch-processing--file-io)
  - [CLI Examples](#cli-examples-1)
  - [Python API Examples](#python-api-examples-1)
- [4. Chemical Structure Downloads](#4-chemical-structure-downloads)
  - [Python API Examples](#python-api-examples-2)
- [5. Advanced SMILES Validation (Python API)](#5-advanced-smiles-validation-python-api)
- [6. Fingerprint Generation](#6-fingerprint-generation)
  - [CLI Examples](#cli)
  - [Python API Examples](#python-api)
- [7. Similarity Search](#7-similarity-search)
  - [CLI Examples](#cli)
  - [Python API Examples](#python-api)
- [8. Substructure Search](#8-substructure-search)\n  - [CLI](#cli)\n  - [Python API](#python-api)\n- [9. Drug-Likeness Filtering (ADMET)](#8-drug-likeness-filtering-admet)
  - [CLI Examples](#cli)
  - [Python API Examples](#python-api)
- [10. Standardization Pipeline](#9-standardization-pipeline)
  - [Python API](#python-api-1)
- [11. Tautomer Enumeration](#10-tautomer-enumeration)
  - [CLI](#cli)
  - [Python API](#python-api-3)
- [12. IUPAC Identifier Generation](#11-iupac-identifier-generation)
  - [Python API](#python-api-2)
- [13. Advanced Cheminformatics](#13-advanced-cheminformatics)
  - [13.1 Reaction SMILES Validation](#131-reaction-smiles-validation)
  - [13.2 Multiple Conformer Generation](#132-multiple-conformer-generation)
  - [13.3 Murcko Scaffold Extraction](#133-murcko-scaffold-extraction)
  - [13.4 Stereochemistry Analysis](#134-stereochemistry-analysis)
  - [13.5 R-Group Decomposition](#135-r-group-decomposition)
  - [13.6 SMILES Augmentation](#136-smiles-augmentation)
  - [13.7 Atom Mapping](#137-atom-mapping)
- [14. SMARTS Pattern Auditing](#14-smarts-pattern-auditing)
  - [14.1 Audit a catalogue](#141-audit-a-catalogue)
  - [14.2 Explain one pattern](#142-explain-one-pattern)
  - [14.3 Declare the preparation when screening](#143-declare-the-preparation-when-screening)
  - [14.4 Python API](#144-python-api)
- [15. Molecular Identity and Library Comparison](#15-molecular-identity-and-library-comparison)
  - [15.1 Identity keys for one molecule](#151-identity-keys-for-one-molecule)
  - [15.2 How many compounds does my file really contain?](#152-how-many-compounds-does-my-file-really-contain)
  - [15.3 Compare two collections](#153-compare-two-collections)
  - [15.4 Python API](#154-python-api)
- [16. SMILES Diagnosis](#16-smiles-diagnosis)
  - [16.1 One string](#161-one-string)
  - [16.2 A file of generated or scraped SMILES](#162-a-file-of-generated-or-scraped-smiles)
  - [16.3 What the checks mean](#163-what-the-checks-mean)
- [Learn More](#learn-more)

---

# 1. System Management

ChemLitmus manages local SQLite caches and logs to ensure high performance and respect for PubChem's rate limits.

### CLI Commands

Check your environment paths, active threads, and database size:

```bash
chemlitmus status
```

Initialize the database schema and storage directories (Run this once after installation):

```bash
chemlitmus init
```

Keep your tool up-to-date. This smart command automatically checks PyPI (or GitHub if you cloned the source) and safely applies updates:

```bash
chemlitmus update
```

If your local database gets corrupted or you want to clear your cache completely, perform a factory reset:

```bash
chemlitmus reinstall -y
```

---

# 2. Single Compound Lookups

ChemLitmus's "Smart Auto-Routing" automatically detects if your input is a SMILES string, an InChIKey, a PubChem CID, or a Chemical Name.

## CLI Examples

### Basic Lookups:

```bash
# Lookup by SMILES
chemlitmus lookup "CC(=O)OC1=CC=CC=C1C(=O)O"

# Lookup by Common/IUPAC Name
chemlitmus lookup "Aspirin"
chemlitmus lookup "benzene"
```

### Advanced CLI Flags:

```bash
# Bypass the local SQLite cache to force a fresh network request
chemlitmus lookup "Caffeine" --no-cache

# Force the engine to treat the input specifically as a CID
chemlitmus lookup 2244 --cid

# Output raw JSON instead of a rich table (ideal for piping into `jq` or other scripts)
chemlitmus lookup "Ibuprofen" --json
```

### Python API Examples

```python
from chemlitmus import lookup, lookup_by_name

# 1. Smart Auto-Detect Lookup
compound = lookup("c1ccccc1")
print(f"Name: {compound.iupac_name}, MW: {compound.molecular_weight}")

# 2. Explicit Lookup by Name (Bypasses SMILES validation checks)
drug = lookup_by_name("Amoxicillin")
print(f"CID: {drug.cid}, Formula: {drug.molecular_formula}")

# 3. Accessing detailed properties (PubChemCompound model)
if drug:
    print(f"H-Bond Donors: {drug.hbond_donor_count}")
    print(f"XLogP: {drug.xlogp}")
    print(f"InChIKey: {drug.inchikey}")
```

---

# 3. Batch Processing & File I/O

Process hundreds of compounds in seconds. ChemLitmus uses multithreading (`concurrent.futures`) combined with a strict rate limiter to fetch data as fast as possible without getting banned by PubChem.

Supported Input Formats: `.csv`, `.tsv`, `.xlsx`, `.smi`, `.sdf`, `.txt`

(ChemLitmus automatically detects the column containing SMILES/Names!)

## CLI Examples

```bash
# Basic batch processing (auto-generates a CSV output and a .log report)
chemlitmus batch input_data.csv

# Output to an Excel file and keep duplicate entries (duplicates are removed by default)
chemlitmus batch raw_smiles.txt --format xlsx --keep-duplicates

# Specify a custom output path and export as JSON
chemlitmus batch data.sdf --output /my_project/clean_data.json --format json
```

## Python API Examples

```python
from chemlitmus import lookup_file, lookup

# 1. Process a file directly in your script
results = lookup_file(
    input_file="messy_data.csv",
    output_file="clean_results.xlsx",
    output_format="xlsx",
    remove_duplicates=True
)

print(f"Successfully processed {len(results)} unique compounds.")

# 2. Custom loop for lists (No file needed)
my_chemicals = ["Aspirin", "c1ccccc1", "Invalid_Chemical_Name"]
valid_compounds = []

for chem in my_chemicals:
    data = lookup(chem, use_cache=True)
    if data and data.cid:
        valid_compounds.append(data)
```

---

# 4. Chemical Structure Downloads & Generation

Download physical structure files from PubChem or generate 2D/3D conformations offline from SMILES using RDKit with built-in resume logic.

Supported Formats: `sdf`, `mol`, `pdb`, `png` (generation supports `sdf`, `mol`, `pdb`)

Supported Dimensions: `2d`, `3d`

Downloaded and generated files are skipped automatically unless `--force` is used.

```bash
# Download a single 3D SDF file by its PubChem CID
chemlitmus download 2244 --format sdf --3d

# Generate all 3D structures locally from SMILES using RDKit (--gen all)
chemlitmus download "CC(=O)OC1=CC=CC=C1C(=O)O" --gen all --3d --format sdf

# Generate 2D MOL structure locally from SMILES
chemlitmus download "c1ccccc1" --gen all --2d --format mol

# Batch download with fallback to local generation (--gen missing)
chemlitmus download -i my_compounds.csv --gen missing --format sdf --3d --output-dir ./structures/

# Batch generate all structures offline from file
chemlitmus download -i my_compounds.smi --gen all --format pdb --3d --output-dir ./3d_models/

# Force overwrite existing files (disables resume logic)
chemlitmus download -i my_compounds.csv --format sdf --force
```

## Python API Examples

```python
from chemlitmus import download_structure

# Download a single structure programmatically
status = download_structure(
    cid=2244, 
    format="sdf", 
    dimension="3d", 
    output_dir="my_structures",
    force=False
)

print(f"Download status: {status}")
```

# 5. Advanced SMILES Validation (Python API)

If you only need to validate SMILES strings and calculate RDKit descriptors locally without querying the PubChem internet database, you can use the core SMILES engine directly.

```python
from chemlitmus import validate_smiles

# Validate a complex SMILES string
result = validate_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")

if result.is_valid:
    print(f"Standardized SMILES: {result.canonical_smiles}")
    print(f"Exact Mass: {result.molecular_weight}")
    print(f"Heavy Atoms: {result.heavy_atom_count}")
    print(f"TPSA: {result.tpsa}")
    print(f"Calculated LogP: {result.logp}")
else:
    print(f"Invalid SMILES! Error: {result.error_message}")
```

---

# Learn More

- **[README.md](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/README.md)** — Installation and Quick Start
- **[API Documentation](https://chemlitmus.readthedocs.io/en/latest/)** — Python API reference
- **[GitHub Issues](https://github.com/AtharvaTilewale/ChemLitmus/issues)** — Bug reports and feature requests
- **[Discussions](https://github.com/AtharvaTilewale/ChemLitmus/discussions)** — Community support

---

# 6. Fingerprint Generation

Generate molecular fingerprints offline from SMILES strings. Useful for ML model preparation, similarity search, and database indexing.

## CLI

```bash
# Single compound - ECFP4 (default)
chemlitmus fingerprint "CC(=O)OC1=CC=CC=C1C(=O)O" --type ecfp4

# MACCS keys (167-bit)
chemlitmus fingerprint "CCO" --type maccs

# All fingerprint types at once
chemlitmus fingerprint "CC(=O)OC1=CC=CC=C1C(=O)O" --type all

# Batch file - save to CSV
chemlitmus fingerprint --file compounds.smi --type ecfp4 --bits 2048 --output fingerprints.csv

# Custom bit size
chemlitmus fingerprint "CC(=O)OC1=CC=CC=C1C(=O)O" --type rdkit --bits 1024
```

**Supported fingerprint types:**

| Type | Algorithm | Default Bits | Use Case |
|------|-----------|-------------|----------|
| `ecfp4` | Morgan (radius=2) | 2048 | General ML, virtual screening |
| `ecfp6` | Morgan (radius=3) | 2048 | More specific substructures |
| `fcfp4` | Feature Morgan (radius=2) | 2048 | Pharmacophore-based |
| `maccs` | MACCS Keys | 167 (fixed) | Structural keys, scaffold analysis |
| `rdkit` | Daylight-style RDKit | 2048 | General purpose |
| `atompair` | Atom Pair | 2048 | 3D-aware searches |
| `torsion` | Topological Torsion | 2048 | Conformer-sensitive searches |

## Python API

```python
from chemlitmus import compute_fingerprint

# Single fingerprint
fp = compute_fingerprint("CC(=O)OC1=CC=CC=C1C(=O)O", fp_type="ecfp4")
print(f"Type: {fp.fingerprint_type}, On bits: {fp.n_on_bits}, Density: {fp.density:.4f}")
print(f"Bit string: {fp.bit_string}")

# All fingerprint types
fps = compute_fingerprint("CCO", fp_type="all")
for fp in fps:
    print(f"{fp.fingerprint_type:12s}: {fp.n_on_bits} bits on / {fp.n_bits}")
```

---

# 7. Similarity Search

Search a compound library to find structurally similar compounds using Tanimoto similarity.

## CLI

```bash
# Search against a library file, default threshold 0.5, top 10
chemlitmus similar "CC(=O)OC1=CC=CC=C1C(=O)O" --file library.smi

# Custom threshold and top-N
chemlitmus similar "CC(=O)OC1=CC=CC=C1C(=O)O" --file library.csv --threshold 0.3 --top 20

# Use MACCS fingerprints instead of ECFP4
chemlitmus similar "CCO" --file compounds.smi --fp-type maccs --top 5

# Save results to CSV
chemlitmus similar "CC(=O)OC1=CC=CC=C1C(=O)O" --file library.smi --output hits.csv

# Higher bit resolution
chemlitmus similar "CCO" --file library.smi --fp-type ecfp6 --bits 4096 --top 10
```

## Python API

```python
from chemlitmus import compute_similarity

# Library as list of SMILES
library = ["CCO", "CCCO", "CC(=O)OC1=CC=CC=C1C(=O)O", "c1ccccc1", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]

hits = compute_similarity(
    query_smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
    library=library,
    fp_type="ecfp4",
    threshold=0.1,
    top_n=5,
)

for hit in hits:
    print(f"Rank {hit.rank}: {hit.hit} (Tanimoto={hit.similarity:.4f})")
```

---


---

# 8. Substructure Search

Search a compound library for molecules that contain a specific substructural fragment or functional group. This is highly useful for identifying compounds with a required pharmacophore or structural alert.

## CLI

```bash
# Search using a SMARTS query (default). This searches for a carboxylic acid.
chemlitmus substructure "C(=O)[OH]" --file library.csv

# Search using a SMILES query (e.g. benzene ring)
chemlitmus substructure "c1ccccc1" --file compounds.smi --smiles-query

# Save only the matching compounds to a new CSV file
chemlitmus substructure "[#9,#17,#35,#53]" --file library.csv --output halogens_only.csv
```

## Python API

```python
from chemlitmus import substructure_search

library = [
    "CCO", 
    "CC(=O)OC1=CC=CC=C1C(=O)O", 
    "c1ccccc1", 
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
]

# Query for carboxylic acid using SMARTS
hits = substructure_search("C(=O)[OH]", library, is_smarts=True)

print(f"Found {len(hits)} matches:")
for hit in hits:
    print(f"- {hit.smiles} (Atoms matched: {hit.match_indices})")
```

# 9. Drug-Likeness Filtering (ADMET)

Evaluate compounds against standard drug-likeness rules and PAINS alerts.

## CLI

```bash
# Single compound - all rules (default)
chemlitmus filter "CC(=O)OC1=CC=CC=C1C(=O)O"

# Specific rules only
chemlitmus filter "CC(=O)OC1=CC=CC=C1C(=O)O" --rules lipinski,veber,pains

# Batch: keep only compounds that pass all rules
chemlitmus filter --file compounds.csv --rules lipinski --output drug_like.csv

# Batch: keep only PAINS-free compounds (--rules pains keeps clean ones)
chemlitmus filter --file compounds.csv --rules pains --output no_pains.csv

# Batch: keep only PAINS-flagged compounds for investigation (--fail inverts)
chemlitmus filter --file compounds.csv --rules pains --fail --output pains_hits.csv

# Add QED minimum threshold
chemlitmus filter --file compounds.csv --rules lipinski,veber --qed-min 0.5 --output output.csv
```

**Available filter rules:**

| Rule | Criteria | Use Case |
|------|----------|----------|
| `lipinski` | MW<=500, LogP<=5, HBD<=5, HBA<=10 | Oral drug candidates |
| `veber` | RotBonds<=10, TPSA<=140 A^2 | Oral bioavailability |
| `ghose` | MW 160-480, LogP -0.4 to 5.6, Atoms 20-70, MR 40-130 | Drug-like space |
| `egan` | TPSA<=131.6, LogP<=5.88 | Passive permeability |
| `ro3` | MW<=300, LogP<=3, HBD<=3, HBA<=3 | Lead-like fragments |
| `pains` | RDKit PAINS catalog | Frequent hitter detection |
| `qed` | 0-1 score (info only) | Overall drug-likeness score |

## Python API

```python
from chemlitmus import apply_filters

# All rules at once
result = apply_filters("CC(=O)OC1=CC=CC=C1C(=O)O")

print(f"MW:  {result.molecular_weight:.2f} g/mol")
print(f"LogP: {result.logp:.2f}")
print(f"QED:  {result.qed_score:.4f}")
print(f"Lipinski: {'PASS' if result.lipinski.passed else 'FAIL'} — {result.lipinski.details}")
print(f"PAINS:    {'PASS' if result.pains.passed else 'FAIL'} — {result.pains.details}")
print(f"Overall:  {'PASS' if result.passes_all else 'FAIL'}")

# Specific rules only
result = apply_filters("CC(=O)OC1=CC=CC=C1C(=O)O", rules=["lipinski", "pains"])

# Batch filtering
import csv

library = ["CC(=O)OC1=CC=CC=C1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "c1ccccc1"]
passing = [smi for smi in library if apply_filters(smi, rules=["lipinski"]).passes_all]
print(f"Drug-like compounds: {len(passing)}/{len(library)}")
```

---

# 10. Standardization Pipeline

Clean up "dirty" SMILES from databases using RDKit's MolStandardize:

```bash
# Single compound: strip salts, neutralize, canonicalize tautomers
chemlitmus standardize "[Na+].[OH-].CC(=O)[O-]"
# Output: CC(=O)O

# Show exactly what each step changed
chemlitmus standardize "[Na+].[OH-].CC(=O)[O-]" --show-diff

# Apply only specific steps
chemlitmus standardize "O=C([O-])c1ccccc1" --steps neutralize,canonical

# Batch process a CSV file
chemlitmus standardize --file compounds.csv --output standardized.csv --show-diff
```

**Available steps** (applied in order):

| Step | Description |
|------|-------------|
| `fragment` | Salt stripping — keeps the largest organic fragment |
| `neutralize` | Neutralizes charged atoms (e.g., carboxylate → carboxylic acid) |
| `tautomer` | Canonicalize tautomers to a single stable form |
| `canonical` | Generate canonical RDKit SMILES string |

### Python API

```python
from chemlitmus import standardize_smiles

result = standardize_smiles("[Na+].[OH-].CC(=O)[O-]")
print(result.output_smiles)   # CC(=O)O
print(result.changed)         # True

# Inspect per-step changes
for step in result.step_results:
    if step.changed:
        print(f"{step.step}: {step.input_smiles} -> {step.output_smiles}")

# Custom pipeline
result = standardize_smiles("O=C([O-])c1ccccc1", steps=["neutralize", "canonical"])
```

---


---

# 11. Tautomer Enumeration

Different tautomers of a single molecule can exhibit vastly different binding affinities to a target protein. Enumerating plausible tautomeric states is a critical preparation step for structure-based virtual screening and molecular docking.

## CLI

```bash
# Single compound
chemlitmus tautomers "Oc1nc(O)c2nc[nH]c2n1"

# Restrict the maximum number of tautomers generated
chemlitmus tautomers "Oc1nc(O)c2nc[nH]c2n1" --max 50

# Batch process a CSV file
# NOTE: The output CSV will "explode" the dataset, meaning if an input SMILES 
# produces 15 tautomers, the output CSV will contain 15 rows for that input.
chemlitmus tautomers --file ligands.csv --output all_tautomers.csv
```

## Python API

```python
from chemlitmus import enumerate_tautomers

# Enumerate tautomers
result = enumerate_tautomers("Oc1nc(O)c2nc[nH]c2n1", max_tautomers=100)

if result.success:
    print(f"Generated {result.num_tautomers} tautomers.")
    print(f"Canonical tautomer: {result.canonical_tautomer}")
    
    # Iterate through them
    for t_smi in result.tautomers:
        is_canonical = (t_smi == result.canonical_tautomer)
        print(f"{t_smi} (Canonical: {is_canonical})")
else:
    print(f"Error: {result.error}")
```

# 12. IUPAC Identifier Generation

Generate systematic identifiers from SMILES, fully offline.

```bash
# Offline: InChI, InChIKey, Formula, MW
chemlitmus iupacname "CC(=O)OC1=CC=CC=C1C(=O)O"

# With IUPAC preferred name (fetched from PubChem on first use, then cached offline)
chemlitmus iupacname "CCO" --online

# Batch processing
chemlitmus iupacname --file compounds.smi --online --output identifiers.csv
```

| Output | Source | Requires Network? |
|--------|--------|-------------------|
| Canonical SMILES | RDKit | No |
| Molecular Formula | RDKit | No |
| Exact MW | RDKit | No |
| InChI | RDKit | No |
| InChIKey | RDKit | No |
| IUPAC Name | PubChem (cached) | First time only |

### Python API

```python
from chemlitmus import get_iupac_name

# Fully offline — InChI, formula, MW
result = get_iupac_name("CCO")
print(result.inchi)        # InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3
print(result.inchikey)     # LFQSCWFLJHTTHZ-UHFFFAOYSA-N
print(result.molecular_formula)  # C2H6O
print(result.molecular_weight)   # 46.0419

# Fetch IUPAC name (network on first call, cache thereafter)
result = get_iupac_name("CCO", use_online=True)
print(result.iupac_name)   # ethanol
print(result.iupac_name_source)  # pubchem (or 'cache' on repeat calls)
```


---

# 13. Advanced Cheminformatics

## 13.1 Reaction SMILES Validation
Validate Reaction SMILES (SMIRKS) syntax, confirming the number of reactants, agents (catalysts), and products.
```bash
chemlitmus reaction "CC(=O)O.OCC>>CC(=O)OCC.O"
```

## 13.2 Multiple Conformer Generation
Generate multiple optimized 3D conformers using the ETKDG algorithm. The output must be saved as a multi-model `.sdf` file, ready for 3D virtual screening.
```bash
chemlitmus conformers "CCO" --num-conformers 100 --output out.sdf
```

## 13.3 Murcko Scaffold Extraction
Remove side chains and extract the core ring framework. Highly useful for clustering High-Throughput Screening (HTS) hits.
```bash
# Single
chemlitmus scaffold "CC(=O)OC1=CC=CC=C1C(=O)O"

# Batch
chemlitmus scaffold --file hits.csv --output scaffolds.csv
```

## 13.4 Stereochemistry Analysis
Identify stereocenters (R/S) and unassigned chiral atoms (`?`). Use `--chiral-flag` to return a system error code if unassigned stereochemistry is detected (useful for strict CI pipelines).
```bash
chemlitmus stereo "C[C@H](O)CC" --chiral-flag
```


## 13.5 R-Group Decomposition
Decompose a library of molecules against a core SMARTS scaffold.
```bash
chemlitmus rgroup --core "c1ccccc1" --smiles "Cc1ccccc1,c1ccccc1F"
```

## 13.6 SMILES Augmentation
Generate a set of unique uncanonical SMILES representing the same molecule, highly useful for deep learning model augmentation.
```bash
chemlitmus augment "CC(=O)OC1=CC=CC=C1C(=O)O" --num 10
```

## 13.7 Atom Mapping
Assign map indices to every atom, useful for reaction tracking and graph networks.
```bash
chemlitmus atommap "CCO"
```


# 14. SMARTS Pattern Auditing

Structural-alert catalogues are redistributed as plain SMARTS text and applied with whatever molecule
preparation a toolkit defaults to. `smartsaudit` evaluates a pattern set against a reference population of
real molecules and reports what is wrong with it before you use it to reject compounds.

## 14.1 Audit a catalogue

```bash
# CSV with a 'smarts' column (optional 'description' and 'rule_set_name'), or a text file of SMARTS
chemlitmus smartsaudit alerts.csv

# Save the per-pattern table and the complete result
chemlitmus smartsaudit alerts.csv --output audit.csv --json audit.json

# Run a subset of checks, tighten the over-broad threshold, list more flagged rows
chemlitmus smartsaudit alerts.csv --checks redundancy,dead --breadth-threshold 0.05 --show 40

# Use your own compound collection as the reference population
chemlitmus smartsaudit alerts.csv --library screening_deck.smi
```

The terminal summary has two tables. The first counts patterns per defect class; the second shows how
many compounds the whole catalogue flags under each molecule preparation and how many compounds change
their pass/fail verdict relative to the default. A catalogue whose verdict flips for a sizeable share of
compounds is not reproducible unless the preparation is declared.

Reading the flags in `audit.csv`:

| Flag | Meaning | Typical action |
|---|---|---|
| `unparseable` | RDKit rejects the SMARTS | Fix the syntax |
| `needs-explicit-h` | Contains a hydrogen atom; silently dead unless molecules carry explicit H | Run with `--prep explicit-h`, or rewrite with H-count primitives (`[CH2]`) |
| `over-broad` | Matches more than the threshold share of the reference set | Check it is a deliberate property filter, not a mis-specified alert |
| `dead:never-matching-atom` | A query atom matches no real atom | Likely defective (e.g. `[N+]#[C-]`, which sanitization never produces) |
| `dead:rare-combination` | Every atom is realisable; the combination does not occur | Probably fine; the alert targets rare chemistry |
| `dead:fires-only-with-<prep>` | Alive under another preparation | Declare that preparation |
| `duplicate` | Identical SMARTS string appears earlier | Remove, or keep for provenance |
| `equivalent` | Identical hit set to another pattern on this reference set | Inspect; may differ on other chemistry |
| `subsumed` | Another pattern's hits contain all of this one's | Redundant as a rejection rule |
| `prep-sensitive` | Hit count changes with preparation | Declare the preparation wherever the pattern is used |

Equivalence and subsumption are empirical over the reference set, not logical properties of the SMARTS.
"Never fires" is an upper bound on dead rules: alert sets deliberately target rare liabilities.

## 14.2 Explain one pattern

```bash
chemlitmus smartsaudit --explain "[$(N(=O)(=O)),$([N+](=O)[O-])]"
chemlitmus smartsaudit --explain "c[H]" --library screening_deck.smi --json explain.json
```

Shows the normalized SMARTS, whether it needs explicit hydrogens, hits under each preparation, example
matches, and for every query atom how many reference molecules contain an atom satisfying that primitive
alone. An atom with zero is the reason a pattern can never fire.

## 14.3 Declare the preparation when screening

```bash
chemlitmus filter --file library.csv --rules pains --prep explicit-h --output screened.csv
chemlitmus substructure "C1=CC=CC=C1" --file library.smi --prep kekule
```

The chosen preparation is written into every output row (`preparation` column), so a screen can be
reproduced exactly.

## 14.4 Python API

```python
from chemlitmus import audit_smarts, explain_smarts, load_patterns, load_reference_library

patterns = load_patterns("alerts.csv")
mols, source = load_reference_library()            # or load_reference_library("screening_deck.smi")

res = audit_smarts(patterns, library=mols, library_source=source, breadth_threshold=0.05)
print(res.n_patterns, res.n_dead, res.n_subsumed, res.n_clean)
print(res.sensitivity.compounds_flagged, res.sensitivity.verdict_flips)

for p in res.patterns:
    if p.dead_verdict == "never-matching atom":
        print(p.name, p.smarts, p.never_matching_atoms)

import csv
rows = res.to_rows()
with open("audit.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

e = explain_smarts("c:1(:c:c:c(:c:c:1)-[#6]=[#7]-[#7])-[#8]-[#1]", library=mols)
print(e.verdict, e.hits_by_preparation)
```


# 15. Molecular Identity and Library Comparison

## 15.1 Identity keys for one molecule

```bash
chemlitmus identity "CC(N)C(=O)O.[Na+].[Cl-]"
```

Shows the key at every level. Levels nest: two records identical at `exact` are identical at every
looser level; two records that first agree at `nostereo` differ only in stereochemistry.

## 15.2 How many compounds does my file really contain?

```bash
chemlitmus identity --file library.csv --level parent
chemlitmus identity --file library.csv --level skeleton --output identity.csv --show 30
```

The first table counts distinct compounds at every level; the second lists groups with more than
one member and what varies among them. The CSV carries every record's keys plus a `group_id` at
the chosen level, so duplicates can be collapsed or inspected downstream.

## 15.3 Compare two collections

```bash
chemlitmus diff chembl_34.smi chembl_35.smi --level parent --output diff.csv
chemlitmus diff vendor_2025.csv vendor_2026.csv --level skeleton --json diff.json
```

Choose the level to match the question. `exact` asks "did any string change"; `parent` ignores
salt and charge form; `nostereo` additionally ignores stereo; `skeleton` ignores stereo and
tautomer form. Compounds present in both collections but written differently are listed as
**changed** with the reason, so a release note can say "41 compounds had stereo assigned, 12 moved
to a different salt form" rather than "1,203 SMILES differ".

## 15.4 Python API

```python
from chemlitmus import compute_identity, group_by_identity, diff_libraries, describe_difference

rep = group_by_identity(smiles, level="parent")
dupes = [g for g in rep.groups if "salt, counter-ion or charge form" in g.differs_by]

d = diff_libraries(old, new, level="parent")
for e in d.entries:
    if e.status == "changed":
        print(e.smiles_a, "->", e.smiles_b, e.change)
```

# 16. SMILES Diagnosis

## 16.1 One string

```bash
chemlitmus diagnose "C1CC(C"
chemlitmus diagnose "CN(C)(C)C"
chemlitmus diagnose "c1cncc1" --no-repair
```

Output: the input with `^` markers under each problem, one line per finding with its category,
message, position or atom index, and a suggestion; then the repaired string and the list of
repairs applied, if a mechanical repair was possible. Exit code 0 for valid, 2 for invalid.

## 16.2 A file of generated or scraped SMILES

```bash
chemlitmus diagnose --file generated.smi --output diagnosis.csv          # invalid records only
chemlitmus diagnose --file generated.smi --output diagnosis.csv --all    # every record
```

The summary table counts invalid records by primary problem and how many were repaired
mechanically. The CSV has one row per record with all problems, positions, suggestions and the
repair outcome.

## 16.3 What the checks mean

| Category | Typical cause |
|---|---|
| characters | a space in the cell (RDKit silently parses only the part before it), an en dash copied from a PDF, a stray letter |
| brackets | `[C@@H` missing its `]`, `[Xx]` unknown element, `[CH3+2-]` malformed |
| parentheses | a branch opened and never closed, or a stray `)` |
| rings | a ring digit opened and never closed (RDKit: "unclosed ring") |
| syntax | anything else RDKit's parser rejects; the position is RDKit's own |
| valence | four-connected neutral nitrogen (needs `[N+]`), five-bonded carbon |
| aromaticity | pyrrole-type `n` without `[nH]`, odd all-carbon aromatic ring, aromatic atom outside a ring |

Repairs are mechanical. They make a string parse; they do not know what molecule was intended.
Treat `repaired_smiles` as a candidate to review, not as the answer.
