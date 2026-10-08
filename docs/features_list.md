# ChemLitmus Features & Commands

ChemLitmus is a production-grade, offline-first command-line tool and Python library for comprehensive cheminformatics workflows.

## Core Features & CLI Commands

### 1. Basic Structure Analysis & Lookup
* **`validate`**: Validate SMILES strings, returning standard properties (Molecular Weight, Formula, InChI, Canonical SMILES).
* **`lookup`**: Query the PubChem REST API by compound name or CID to fetch structures and identifiers.
* **`download`**: Generate 2D or 3D structures (SDF, MOL, PDB, XYZ) for a single SMILES or a batch file.

### 2. Advanced Cheminformatics
* **`reaction`**: Parse, validate, and analyze Reaction SMILES (SMIRKS). Checks atom mapping, reactants, agents (catalysts), and products.
* **`conformers`**: Generate multiple optimized 3D conformer ensembles (ETKDG + MMFF94) for virtual screening, outputting to multi-model `.sdf`.
* **`scaffold`**: Extract the Murcko Scaffold framework from molecules, stripping side-chains. Supports batch processing for HTS library clustering.
* **`stereo`**: Analyze stereocenters (R/S) and flag unassigned chiral centers (`?`). Includes strict CI/CD integration via `--chiral-flag`.

### 3. Molecular Standardization
* **`standardize`**: Clean up "dirty" SMILES through salt stripping, neutralization, tautomer canonicalization, and fragment removal.
* **`iupacname`**: Fetch and generate IUPAC names using a hybrid offline InChI-based algorithm and a local SQLite cache.
* **`tautomers`**: Enumerate all valid tautomeric states of a molecule, with batch processing for AutoDock Vina preparation.

### 4. Search, Filtering, and Similarity
* **`fingerprint`**: Compute bit-vector molecular fingerprints (MACCS, ECFP4, ECFP6) for machine learning and similarity scoring.
* **`similar`**: Search a local library (CSV/SMI) to find molecules similar to a query SMILES based on Tanimoto similarity thresholds.
* **`substructure`**: Perform strict Substructure Searches against a library using SMARTS patterns or exact SMILES fragments.
* **`filter`**: Apply strict ADMET drug-likeness rules (e.g., Lipinski's Rule of 5) to screen out undesirable compounds from a dataset.

### 5. Utilities
* **`update`**: Automatically download and install the latest version of ChemLitmus from GitHub or PyPI.

---
## Python API Features
All CLI commands are fully exposed as a typed Python API in the `chemlitmus.core` namespace, returning `pydantic` models (e.g., `SMILESValidationResult`, `FilterResult`, `ScaffoldResult`) for robust programmatic integration into data pipelines.
