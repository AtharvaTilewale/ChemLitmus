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

- **`smartsproof`** — static SMARTS containment proofs with no reference library: atom and bond expressions are evaluated over a finite universe of atom states built from the patterns themselves, and `P ⊆ Q` is proven by an injective witness mapping (Chandra–Merlin criterion under RDKit's match semantics). Proves redundancy within a catalogue, decides one pair (`--subsumes`, `--equivalent`), and checks local satisfiability (`--satisfiable`). Sound, not complete: *no witness* is reported as undecided, never as refuted; recursive SMARTS, isotopes, valence and ring-size primitives are reported as not analysable. `smartsaudit --checks all,proof` adds `proven_subsumed_by` / `proof_status` beside the empirical columns.
- **`smartsdiff`** — semantic diff of two SMARTS catalogues: patterns paired by name, identical SMARTS or identical hit set; each pair classified as same hits / broadened / narrowed / shifted / broken / repaired, unpaired as added / removed; the number of reference molecules whose flagged verdict changes between the catalogues. Either side may be an RDKit built-in catalogue (`rdkit:PAINS`, `rdkit:BRENK`, …), compared on hits alone.

- **`smartsaudit`** — audit a SMARTS pattern set (structural alerts, substructure filters) against a reference molecule population: unparseable patterns, patterns that need explicit hydrogens, over-broad patterns, dead patterns triaged into never-matching atom / rare combination / fires-only-under-another-preparation, exact and library-equivalent duplicates, strict subsumption, and sensitivity of hit counts and per-compound verdicts to molecule preparation (implicit H, explicit H, kekulized). `--explain` decomposes a single pattern atom by atom. Matching runs multithreaded through RDKit's `SubstructLibrary`.
- Bundled 9,272-molecule ChEMBL-derived reference library (`chemlitmus/data/`, CC BY-SA 3.0); `--library` substitutes any SMILES-bearing file.
- **`--prep`** on `filter` and `substructure` — declare the molecule preparation used for substructure matching; recorded in every output row (`FilterResult.preparation`, `SubstructureHit.preparation`).
- Python API: `audit_smarts`, `explain_smarts`, `load_patterns`, `load_reference_library`, `prepare_molecule`; models `SmartsAuditResult`, `PatternAudit`, `SensitivitySummary`, `SmartsExplanation`, `AtomExplanation`.

### Validation

- **`validate`** — offline SMILES validation with canonical form, formula, exact mass, LogP, HBD/HBA, TPSA and heavy-atom count; `--quiet` prints only the canonical SMILES; exit code 2 on an invalid single SMILES. Whitespace inside a SMILES is reported as invalid, because RDKit would otherwise silently parse only the first token.
- **`init`** — create the cache, data and log directories and the SQLite database, and report their locations.

### Data access and infrastructure

- **`resolve`** — one query across **PubChem, ChEMBL, ChEBI and KEGG** (opt-in) in parallel under a single time budget; every database answers in one `CompoundRecord` schema. Reports whether the sources **agree on the structure** (InChIKey majority), merges the records field-wise, and collects cross-database identifiers (CAS, DrugBank, HMDB, PDBe, SureChEMBL, …) through **UniChem** and the databases' own links. A source that disagrees with the majority, or does not know the name, is re-queried by the consensus InChIKey, so a free-text search that lands on a derivative is corrected by structure.
- **`concordance`** — for a list of names, measure how often the databases return the same structure, at which identity level they part ways (salt form, tautomer, stereochemistry, different compound), and how often each database's text search landed on a different structure than the consensus.
- `resolve` explains disagreements in identity-level language (`agreement_level`, `disagreement`, `pairwise`) and reports every structure-based correction of a text hit; the consensus on a tie prefers the source that knows the name verbatim, then the parent form over a salt or hydrate. KEGG records now carry SMILES/InChI/InChIKey computed from the MOL block, so KEGG takes part in the agreement check; KEGG name search requires an exact name (its `find` is a substring search) and falls back to KEGG DRUG. ChEBI names are stripped of presentational HTML.
- Local SQLite cache for every database provider: a record found once is reusable offline by query, native identifier or InChIKey (`resolve --no-cache` bypasses it).
- Provider layer `chemlitmus.providers` — `resolve()`, `get_provider()`, `unichem_xrefs()`, `PubChemProvider`, `ChEMBLProvider`, `ChEBIProvider`, `KEGGProvider`; per-database rate limiting and deadline-aware retries so a busy server cannot exceed the fan-out budget.
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
- Input parsing: a header-less list of names or identifiers (e.g. `aspirin`, `CHEMBL25`, `2244`) no longer loses its first line — the first row is treated as a header only when it is a recognised column name; files with a UTF-8 BOM are read correctly.
- Removed unused imports; `ruff` configuration migrated to the `[tool.ruff.lint]` table.

### Tests

- Offline, HTTP-mocked tests for the PubChem client, result export, UniChem and the provider cache; every test runs against an isolated temporary cache.
- CI runs lint, tests on Python 3.10–3.13 (Linux, macOS, Windows), documentation build and package build; live-database tests run on push and are non-blocking.

- 215 offline tests plus 8 networked integration tests, running on Python 3.10–3.12.
