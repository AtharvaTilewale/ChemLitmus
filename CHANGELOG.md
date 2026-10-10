# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-10-09

Initial release.

### Dataset quality and evaluation

- **`audit`** — one offline pass over a dataset: parse status, policy-driven standardisation with before/after provenance (text change vs *identity* change), component / charge / isotope / stereo / element flags, descriptors, configurable structural alerts with matched atoms, identity keys at all six levels and groups at the policy level, split leakage and endpoint conflicts. 32 stable issue codes, each with severity, affected records, evidence and a suggested action. Outputs: annotated records CSV, issues CSV, identity groups CSV, summary JSON, authoritative `audit.json`, resolved policy, run manifest, an offline HTML report (escaped, structure depictions with alert highlights, descriptor histograms, searchable issue table), and an optional clean export whose every exclusion is listed. Input accounting reconciles exactly.
- **`leakage`** — overlap between splits as distinct, non-additive evidence classes (exact, parent, tautomer, nostereo, skeleton), with scaffold overlap and nearest-neighbour similarity reported separately as *relatedness*, formula matches reported and explicitly not called identity, within-split duplicates, and an optional temporal-order check.
- **`conflicts`** — contradictory labels or measurements among records of the same compound within the same endpoint context; molar unit conversion with provenance; censored values (`<`, `>`, ranges) kept as bounds and never averaged; tolerance in log10 units; default action is review.
- **`cliffs`** — activity cliffs and label outliers among near-identical *different* compounds: matched molecular pairs (single cut plus hydrogen positions; varied fragment ≤ 13 heavy atoms and never larger than the core) and fingerprint neighbours at Tanimoto ≥ 0.9, compared within endpoint context after same-compound records are merged; a pair is a cliff only when the measurements prove a difference at or above the threshold, censored values make it undetermined; per-transformation summaries with signed mean differences; label outliers — compounds whose every neighbour is a proven cliff against them while at least two of those neighbours agree with each other (`n_neighbours` / `n_agreeing` reported) — listed for review, never corrected. Issue codes `ACTIVITY_CLIFF`, `LABEL_OUTLIER`.
- **`split`** — reproducible random / identity / scaffold / temporal / source splits that never divide a group; acyclic molecules each form their own scaffold group; reports achieved fractions, group-size distribution, endpoint balance and every unsatisfiable constraint; verified with the same leakage implementation.
- **`generated`** — evaluation of generated molecules with explicit denominators: validity over all attempts, uniqueness over valid outputs at each identity level, novelty over unique valid outputs against the named reference sets only; scaffold diversity, nearest-reference similarity, alert prevalence, explicit constraints; repaired candidates evaluated as a separate population; metrics without a denominator reported as `undefined` rather than zero.
- **Record-preserving ingestion** (`chemlitmus.core.records`) — CSV/TSV/XLSX/SMI/SDF into records with stable ids, positions, source ids, original fields, status and issues; ambiguous column mapping is a `SchemaError`, never a guess; SDF properties and failed records preserved.
- **Versioned chemical policy** (`ChemicalPolicy`) — validated and hashed; presets `conservative` and `parent`; controls fragments, charges, tautomer/stereo/isotope handling, identity level, match preparation, alert sets, fingerprint, repair mode and per-code severities.
- **Run manifests** (`RunManifest`) — input/output checksums, tool / Python / RDKit versions, policy hash and settings, written with every audit.
- **Gates** — `--fail-on-severity`, `--max-invalid-fraction`, `--max-duplicate-fraction`, `--no-split-overlap`, `--no-label-conflicts`, `--require-complete`; exit 0 pass, 1 configuration error, 3 policy violation, 4 partial processing.
- **Technical vs biological variability** — with a source/document column, each conflict group reports how much of its disagreement lies within one source (repeated entries from one experiment) and how much lies between sources; without that metadata the split is reported as `not determinable` rather than guessed.
- **Leakage scope** — the report states that it measures contamination between the collections supplied to it and cannot see a model's pretraining corpus, so an absence of overlap is not evidence that a model never saw these compounds.
- **Provider retrieval provenance** — `CompoundRecord` carries `retrieved_at`, `from_cache`, `cached_at` and (on merged records) `field_provenance` naming the source of every merged value, so stale cache data is distinguishable from a fresh fetch.
- **Parquet input** as an optional extra (`pip install 'chemlitmus[parquet]'`).
- Supported Python (3.10–3.13) and RDKit (>=2023.09) ranges documented, with the reason RDKit versions matter.

### Validation and diagnosis

- **`validate`** — offline SMILES validation with canonical form, formula, exact mass, LogP, HBD/HBA, TPSA and heavy-atom count; `--quiet` prints only the canonical SMILES; exit code 2 on an invalid single SMILES. Whitespace inside a SMILES is reported as invalid, because RDKit would otherwise silently parse only the first token.
- **`init`** — create the cache, data and log directories and the SQLite database, and report their locations.

### Standardisation, identity and comparison

- **`validate`** — SMILES validation and canonicalization via RDKit, returning molecular weight, formula, LogP, HBD/HBA, TPSA, heavy-atom count and InChI.
- **`standardize`** — offline standardization pipeline built on RDKit `MolStandardize`: salt stripping (largest-fragment selection), charge neutralization, tautomer canonicalization and canonical SMILES output. Steps are individually selectable, and `--show-diff` reports what changed at each stage.
- **`tautomers`** — enumeration of plausible tautomeric states, with batch processing that emits one row per tautomer.
- **`stereo`** — stereocenter analysis reporting R/S assignments and flagging unassigned centers, with `--chiral-flag` for CI enforcement.
- **`iupacname`** — offline InChI, InChIKey, molecular formula and exact mass; optional IUPAC systematic name from PubChem with local SQLite caching.
- **`atommap`** — assignment of unique atom map numbers.
- **`identity`** — layered molecular identity keys (`exact` > `parent` > `tautomer` / `nostereo` > `skeleton` > `formula`, built on RDKit `RegistrationHash`). Single-molecule key table, or group a collection at any level with per-level distinct-compound counts, multi-member groups and what varies within each (salt/charge form, tautomer, stereochemistry). CSV export with keys and group ids.
- **`diff`** — structure-aware comparison of two compound collections at a chosen identity level: added, removed, unchanged, and changed-with-reason, plus multiplicity changes and Jaccard overlap. CSV and JSON export.
- **`diagnose`** — deterministic, explainable SMILES failure diagnosis. Seven ordered checks (characters, bracket atoms, parentheses, ring closures, RDKit syntax, valence, aromaticity), each located to a character position or atom with a suggestion; safe mechanical repairs (whitespace and dash normalisation, dangling branches, unclosed ring digits, `[nH]`, non-ring aromatic atoms) are attempted and re-validated. Detects records that RDKit silently truncates at whitespace. Exit code 2 on an invalid single SMILES.
- Python API: `compute_identity`, `group_by_identity`, `strictest_shared_level`, `describe_difference`, `diff_libraries`, `diagnose_smiles`; models `IdentityKeys`, `IdentityGroup`, `IdentityReport`, `DiffEntry`, `LibraryDiff`, `SmilesProblem`, `SmilesDiagnosis`; constants `IDENTITY_LEVELS`, `DIAGNOSTIC_CATEGORIES`.

### Pattern quality control

- **Catalogue-cleanup proposals** (`smartsaudit --cleanup`) — each proposal carries its evidence and one of three actions: `recommended` (proof-backed, removal cannot change a verdict), `review` (observed on this panel only) or `keep provenance` (the only cover is in a different published set). Nothing is removed automatically.
- **`smartsproof`** — static SMARTS containment proofs with no reference library: atom and bond expressions are evaluated over a finite universe of atom states built from the patterns themselves, and `P ⊆ Q` is proven by an injective witness mapping (Chandra–Merlin criterion under RDKit's match semantics). Proves redundancy within a catalogue, decides one pair (`--subsumes`, `--equivalent`), and checks local satisfiability (`--satisfiable`). Sound, not complete: *no witness* is reported as undecided, never as refuted; recursive SMARTS, isotopes, valence and ring-size primitives are reported as not analysable. `smartsaudit --checks all,proof` adds `proven_subsumed_by` / `proof_status` beside the empirical columns.
- **`smartsdiff`** — semantic diff of two SMARTS catalogues: patterns paired by name, identical SMARTS or identical hit set; each pair classified as same hits / broadened / narrowed / shifted / broken / repaired, unpaired as added / removed; the number of reference molecules whose flagged verdict changes between the catalogues. Either side may be an RDKit built-in catalogue (`rdkit:PAINS`, `rdkit:BRENK`, …), compared on hits alone.

- **`smartsaudit`** — audit a SMARTS pattern set (structural alerts, substructure filters) against a reference molecule population: unparseable patterns, patterns that need explicit hydrogens, over-broad patterns, dead patterns triaged into never-matching atom / rare combination / fires-only-under-another-preparation, exact and library-equivalent duplicates, strict subsumption, and sensitivity of hit counts and per-compound verdicts to molecule preparation (implicit H, explicit H, kekulized). `--explain` decomposes a single pattern atom by atom. Matching runs multithreaded through RDKit's `SubstructLibrary`.
- Bundled 9,272-molecule ChEMBL-derived reference library (`chemlitmus/data/`, CC BY-SA 3.0); `--library` substitutes any SMILES-bearing file.
- **`--prep`** on `filter` and `substructure` — declare the molecule preparation used for substructure matching; recorded in every output row (`FilterResult.preparation`, `SubstructureHit.preparation`).
- Python API: `audit_smarts`, `explain_smarts`, `load_patterns`, `load_reference_library`, `prepare_molecule`; models `SmartsAuditResult`, `PatternAudit`, `SensitivitySummary`, `SmartsExplanation`, `AtomExplanation`.
- `smartsaudit` records **evidence labels** distinguishing exact text duplicates, observed hit-set identity/containment (claims about *this panel*), static proven containment/equivalence (universal), "not observed in this reference", "unsupported" and "undecided"; the **match semantics** (preparation, hydrogen handling, chirality ignored as RDKit does by default, RDKit version); **catalogue metadata** (`--catalogue-name/-version/-source/-licence`) with a content hash and the pattern count per rule set; a **reference-panel summary** (size, median heavy atoms, ring and charge fractions, element census); **example matches with matched atoms** and preparation-flip examples; and **`--holdout`** libraries, which name the patterns that never fire on the main panel but do fire elsewhere.
- `benchmarks/alert_preparation_sensitivity.py` — a versioned recipe that reproduces the preparation-sensitivity figures per rule set with the catalogue checksum.

### Screening and search

- **`fingerprint`** — ECFP4, ECFP6, FCFP4, MACCS, RDKit, AtomPair and topological torsion fingerprints, offline.
- **`similar`** — Tanimoto similarity search against a local library with configurable threshold and top-N ranking.
- **`substructure`** — substructure search using SMARTS patterns or strict SMILES queries.
- **`filter`** — drug-likeness and ADMET screening: Lipinski Ro5, Veber, Ghose, Egan, Rule of Three, PAINS alerts and QED scoring, with `--fail` to invert selection.
- **`scaffold`** — Murcko scaffold extraction, batch-capable for library clustering.
- **`rgroup`** — R-group decomposition against a common core SMARTS.

### Structures and 3D

- **`download`** — 2D/3D structure retrieval from PubChem in SDF, MOL, PDB and PNG, with resume logic.
- **`--gen`** — offline 2D/3D structure generation from SMILES via RDKit with MMFF94/UFF optimization; `--gen missing` falls back to local generation when PubChem has no record.
- **`conformers`** — multi-conformer ensemble generation (ETKDG + MMFF94) written as multi-model SDF.
- **`reaction`** — reaction SMILES/SMIRKS parsing and validation, reporting reactant, agent and product counts.
- **`augment`** — randomized non-canonical SMILES generation for machine-learning data augmentation.

### Data access

- **`resolve`** — one query across **PubChem, ChEMBL, ChEBI and KEGG** (opt-in) in parallel under a single time budget; every database answers in one `CompoundRecord` schema. Reports whether the sources **agree on the structure** (InChIKey majority), merges the records field-wise, and collects cross-database identifiers (CAS, DrugBank, HMDB, PDBe, SureChEMBL, …) through **UniChem** and the databases' own links. A source that disagrees with the majority, or does not know the name, is re-queried by the consensus InChIKey, so a free-text search that lands on a derivative is corrected by structure.
- **`concordance`** — for a list of names, measure how often the databases return the same structure, at which identity level they part ways (salt form, tautomer, stereochemistry, different compound), and how often each database's text search landed on a different structure than the consensus.
- Disagreements are explained in identity-level language (`agreement_level`, `disagreement`, `pairwise`) and reports every structure-based correction of a text hit; the consensus on a tie prefers the source that knows the name verbatim, then the parent form over a salt or hydrate. KEGG records carry SMILES/InChI/InChIKey computed from the MOL block, so KEGG takes part in the agreement check; KEGG name search requires an exact name (its `find` is a substring search) and falls back to KEGG DRUG. ChEBI names are stripped of presentational HTML.
- Local SQLite cache for every database provider: a record found once is reusable offline by query, native identifier or InChIKey (`resolve --no-cache` bypasses it).
- Provider layer `chemlitmus.providers` — `resolve()`, `get_provider()`, `unichem_xrefs()`, `PubChemProvider`, `ChEMBLProvider`, `ChEBIProvider`, `KEGGProvider`; per-database rate limiting and deadline-aware retries so a busy server cannot exceed the fan-out budget.
- **`lookup`** / **`batch`** — PubChem queries by SMILES, CID, name, InChI or InChIKey, with automatic routing and name-based fallback.
- Multi-format input parsing for CSV, TSV, XLSX, SMI, SDF and TXT with SMILES column auto-detection.
- Multithreaded batch processing with thread-safe rate limiting and retry logic.
- Local SQLite caching of PubChem results.
- Export to CSV, Excel and JSON.
- Typed Python API: every core function returns a Pydantic model.
- **`status`**, **`init`** and **`update`** utility commands.

### Semantics worth knowing

- **Fingerprint names describe the computation.** `fcfp4` uses feature-class (pharmacophoric) Morgan atom invariants and `ecfp4` connectivity invariants, verified against RDKit's own generators; `FingerprintResult` records `algorithm`, `radius`, `atom_invariants`, `representation` and `chirality`.
- **One strict structure parser.** `mol_from_smiles()` backs every structure-consuming operation, so a field with internal whitespace is an error rather than the first token: `CC O` is never silently read as ethane. CXSMILES is parsed as such, and `.smi` lines — which *define* a name after whitespace — are split explicitly.
- **A failed computation is never a chemical negative.** `smartsaudit` records `match_error` per pattern (hit counts unknown, not zero) and a `PreparationStatus` per preparation (a molecule that cannot be Kekulised is *unevaluated*, not a non-hit, and is excluded from verdict-flip counts); `smartsdiff` excludes unprepared molecules and reports how many; the dataset audit reports `processing_complete`.
- **Repairs are candidates, not results.** `diagnose` returns `repair_status` (`not attempted` / `not needed` / `candidate` / `failed`). A candidate parses; it is not claimed to be the intended molecule, carries no confidence number, and is applied only under an explicit repair policy, with the original record always retained.
- **Identity levels nest and their counts are not additive.** `tautomer` and `nostereo` are siblings below `parent`, so the difference between consecutive rows of an identity report is not a count of one cause; `differs_by` names the axis per group. Equal formulas identify a formula group, not a compound.
- **Endpoints without a units column are dimensionless, not incomparable.** `conflicts` and `cliffs` compare a pIC50, logP or score column on its own scale when no units column is supplied; when one is supplied, a value lacking units is flagged `UNITS_MISSING` and excluded rather than guessed.
- **Alert measurements name their denominator and their catalogue.** Preparation-sensitivity figures are reported per rule set, never as a union figure attributed to one named set: on the bundled 9,272-molecule panel PAINS flags 3.8% / 5.0% / 2.1% of compounds under implicit-H / explicit-H / Kekulé (1.3% verdict flips) while SureChEMBL flags 26.2% / 68.7% / 41.2% (42.5% flips); the 77.6% / 90.5% / 88.5% figures cover all eight published sets together. `benchmarks/alert_preparation_sensitivity.py` reproduces every figure with the catalogue checksum.

### Command behaviour

- `conflicts --source-column` separates within-source (technical) from between-source spread per conflict group, as `audit --source-column` already did; the terminal table shows both spreads and the number of sources.
- Invalid option values (`--prep`, `--rules`, `--steps`, `--type`, `--checks`, `--level`) are usage errors (exit 1) rejected before any work or network request, never counted as per-compound failures.
- `conformers` accepts `--num` as an alias for `--num-conformers` and defaults `--output` to `conformers.sdf`; `rgroup` takes the core SMARTS positionally or with `--core`; `scaffold` reports acyclic molecules as `acyclic`.
- `download --gen all` is fully offline: it never contacts PubChem for a title lookup.
- `diagnose` routes RDKit parse messages into the Python logger, so captured diagnostics do not leak to stderr for later parses.
- Input parsing treats the first row as a header only when a cell is a recognised column name, so a header-less list of names or identifiers (`aspirin`, `CHEMBL25`, `2244`) keeps its first line; UTF-8 BOMs are handled.

### Tests and tooling

- **Example data** — `examples/` holds small hand-built inputs with planted, documented defects (`bioactivity.csv`, `library.csv` / `library_v2.csv` / `library.smi`, `alerts.csv` / `alerts_v2.csv`, `generated.smi` / `train.smi`, `demo.csv`, `policy.json`, …); every input file named in the documentation is one of them, `examples/README.md` says what each command finds, and `tests/test_examples.py` runs the documented commands against them so the docs cannot drift from the tool.
- 376 offline tests and 10 networked integration tests. Offline tests mock HTTP at the provider boundary and run against an isolated temporary cache, so the suite never touches the user's cache or the network.
- A curated regression corpus covering invalid, truncated and whitespace-bearing representations, quoted CSV fields, blank rows, failed SDF records, salts, mixtures, zwitterions, permanent charges, isotopes, tautomers, assigned and partial stereochemistry, acyclic molecules and out-of-scope chemistry — each case with the reason it is expected to behave as it does.
- Invariant checks: standardisation idempotence, identity invariance under equivalent SMILES spellings, group-preserving splits, exact input accounting, and CLI/API parity on the same input.
- CI runs lint, the offline suite on Python 3.10–3.13 (Linux, plus macOS and Windows on 3.12), a documentation build and a package build; live-database tests run on push and are non-blocking.
