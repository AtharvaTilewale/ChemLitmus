# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-10-08

Initial release.

### Structure handling

- **`validate`** — SMILES validation and canonicalization via RDKit, returning molecular weight, formula, LogP, HBD/HBA, TPSA, heavy-atom count and InChI.
- **`standardize`** — offline standardization pipeline built on RDKit `MolStandardize`: salt stripping (largest-fragment selection), charge neutralization, tautomer canonicalization and canonical SMILES output. Steps are individually selectable, and `--show-diff` reports what changed at each stage.
- **`tautomers`** — enumeration of plausible tautomeric states, with batch processing that emits one row per tautomer.
- **`stereo`** — stereocenter analysis reporting R/S assignments and flagging unassigned centers, with `--chiral-flag` for CI enforcement.
- **`iupacname`** — offline InChI, InChIKey, molecular formula and exact mass; optional IUPAC systematic name from PubChem with local SQLite caching.
- **`atommap`** — assignment of unique atom map numbers.

### Search and filtering

- **`fingerprint`** — ECFP4, ECFP6, FCFP4, MACCS, RDKit, AtomPair and topological torsion fingerprints, offline.
- **`similar`** — Tanimoto similarity search against a local library with configurable threshold and top-N ranking.
- **`substructure`** — substructure search using SMARTS patterns or strict SMILES queries.
- **`filter`** — drug-likeness and ADMET screening: Lipinski Ro5, Veber, Ghose, Egan, Rule of Three, PAINS alerts and QED scoring, with `--fail` to invert selection.
- **`scaffold`** — Murcko scaffold extraction, batch-capable for library clustering.
- **`rgroup`** — R-group decomposition against a common core SMARTS.

### Structure generation

- **`download`** — 2D/3D structure retrieval from PubChem in SDF, MOL, PDB and PNG, with resume logic.
- **`--gen`** — offline 2D/3D structure generation from SMILES via RDKit with MMFF94/UFF optimization; `--gen missing` falls back to local generation when PubChem has no record.
- **`conformers`** — multi-conformer ensemble generation (ETKDG + MMFF94) written as multi-model SDF.
- **`reaction`** — reaction SMILES/SMIRKS parsing and validation, reporting reactant, agent and product counts.
- **`augment`** — randomized non-canonical SMILES generation for machine-learning data augmentation.

### Data access and infrastructure

- **`lookup`** / **`batch`** — PubChem queries by SMILES, CID, name, InChI or InChIKey, with automatic routing and name-based fallback.
- Multi-format input parsing for CSV, TSV, XLSX, SMI, SDF and TXT with SMILES column auto-detection.
- Multithreaded batch processing with thread-safe rate limiting and retry logic.
- Local SQLite caching of PubChem results.
- Export to CSV, Excel and JSON.
- Typed Python API: every core function returns a Pydantic model.
- **`status`**, **`init`** and **`update`** utility commands.

### Tests

- 136 offline tests plus 8 networked integration tests, running on Python 3.10–3.12.
