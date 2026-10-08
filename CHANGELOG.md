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

### Identity, comparison and diagnosis

- **`identity`** — layered molecular identity keys (`exact` > `parent` > `tautomer` / `nostereo` > `skeleton` > `formula`, built on RDKit `RegistrationHash`). Single-molecule key table, or group a collection at any level with per-level distinct-compound counts, multi-member groups and what varies within each (salt/charge form, tautomer, stereochemistry). CSV export with keys and group ids.
- **`diff`** — structure-aware comparison of two compound collections at a chosen identity level: added, removed, unchanged, and changed-with-reason, plus multiplicity changes and Jaccard overlap. CSV and JSON export.
- **`diagnose`** — deterministic, explainable SMILES failure diagnosis. Seven ordered checks (characters, bracket atoms, parentheses, ring closures, RDKit syntax, valence, aromaticity), each located to a character position or atom with a suggestion; safe mechanical repairs (whitespace and dash normalisation, dangling branches, unclosed ring digits, `[nH]`, non-ring aromatic atoms) are attempted and re-validated. Detects records that RDKit silently truncates at whitespace. Exit code 2 on an invalid single SMILES.
- Python API: `compute_identity`, `group_by_identity`, `strictest_shared_level`, `describe_difference`, `diff_libraries`, `diagnose_smiles`; models `IdentityKeys`, `IdentityGroup`, `IdentityReport`, `DiffEntry`, `LibraryDiff`, `SmilesProblem`, `SmilesDiagnosis`; constants `IDENTITY_LEVELS`, `DIAGNOSTIC_CATEGORIES`.

### Pattern quality control

- **`smartsaudit`** — audit a SMARTS pattern set (structural alerts, substructure filters) against a reference molecule population: unparseable patterns, patterns that need explicit hydrogens, over-broad patterns, dead patterns triaged into never-matching atom / rare combination / fires-only-under-another-preparation, exact and library-equivalent duplicates, strict subsumption, and sensitivity of hit counts and per-compound verdicts to molecule preparation (implicit H, explicit H, kekulized). `--explain` decomposes a single pattern atom by atom. Matching runs multithreaded through RDKit's `SubstructLibrary`.
- Bundled 9,272-molecule ChEMBL-derived reference library (`chemlitmus/data/`, CC BY-SA 3.0); `--library` substitutes any SMILES-bearing file.
- **`--prep`** on `filter` and `substructure` — declare the molecule preparation used for substructure matching; recorded in every output row (`FilterResult.preparation`, `SubstructureHit.preparation`).
- Python API: `audit_smarts`, `explain_smarts`, `load_patterns`, `load_reference_library`, `prepare_molecule`; models `SmartsAuditResult`, `PatternAudit`, `SensitivitySummary`, `SmartsExplanation`, `AtomExplanation`.

### Validation

- **`validate`** — offline SMILES validation with canonical form, formula, exact mass, LogP, HBD/HBA, TPSA and heavy-atom count; `--quiet` prints only the canonical SMILES; exit code 2 on an invalid single SMILES. Whitespace inside a SMILES is reported as invalid, because RDKit would otherwise silently parse only the first token.
- **`init`** — create the cache, data and log directories and the SQLite database, and report their locations.

### Data access and infrastructure

- **`lookup`** / **`batch`** — PubChem queries by SMILES, CID, name, InChI or InChIKey, with automatic routing and name-based fallback.
- Multi-format input parsing for CSV, TSV, XLSX, SMI, SDF and TXT with SMILES column auto-detection.
- Multithreaded batch processing with thread-safe rate limiting and retry logic.
- Local SQLite caching of PubChem results.
- Export to CSV, Excel and JSON.
- Typed Python API: every core function returns a Pydantic model.
- **`status`**, **`init`** and **`update`** utility commands.

### Fixed

- `filter`: an invalid `--prep` or `--rules` value is now a usage error (exit 1) instead of being counted silently as a per-compound failure.
- `standardize`: an invalid `--steps` value no longer prints a second, empty `Error:` line.
- `tautomers`: exits 1 on an invalid single SMILES instead of 0.
- `lookup`: an unknown `--type` is rejected before any network request.
- `conformers`: `--num` accepted as an alias for `--num-conformers`; `--output` defaults to `conformers.sdf`.
- `rgroup`: the core SMARTS may be given positionally as well as with `--core`.
- `scaffold`: acyclic molecules are reported as `acyclic` instead of an empty cell.
- `download --gen all`: no longer contacts PubChem for a title lookup; the command is fully offline as documented.
- `diagnose`: RDKit parse messages captured during diagnosis no longer leak to stderr for subsequent parses.
- Removed unused imports; `ruff` configuration migrated to the `[tool.ruff.lint]` table.

### Tests

- 215 offline tests plus 8 networked integration tests, running on Python 3.10–3.12.
