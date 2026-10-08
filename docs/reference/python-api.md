# Python API reference

Everything is importable from the top-level package:

```python
import chemlitmus
from chemlitmus import validate_smiles, diagnose_smiles, audit_smarts  # etc.
```

Conventions:

- Every function returns a **Pydantic model** (or a list of them). Models are immutable-by-convention records with `.model_dump()` / `.model_dump_json()` for serialisation and `Model.model_validate(dict)` to rebuild.
- **Errors are values.** A bad SMILES does not raise; the result carries `error` (or `error_message`, `is_valid=False`). Only programming errors — an unknown level name, an unknown preparation — raise `ValueError`.
- Functions that take a `preparation` or `level` accept the strings listed in the constants below.

## Constants

| Name | Value |
|---|---|
| `IDENTITY_LEVELS` | `['exact', 'parent', 'tautomer', 'nostereo', 'skeleton', 'formula']` |
| `PREPARATIONS` | `['implicit-h', 'explicit-h', 'kekule']` |
| `AUDIT_CHECKS` | `['compile', 'breadth', 'dead', 'redundancy', 'sensitivity']` |
| `DIAGNOSTIC_CATEGORIES` | `['characters', 'brackets', 'parentheses', 'rings', 'syntax', 'valence', 'aromaticity']` |
| `STANDARDIZE_STEPS` | `['fragment', 'neutralize', 'tautomer', 'canonical']` |

---

## Validation and diagnosis

### `validate_smiles(smiles_str) -> SMILESValidationResult`

Offline validation with basic descriptors. Whitespace inside the string is reported as invalid.

`SMILESValidationResult`: `input_smiles, is_valid, canonical_smiles, molecular_formula, molecular_weight` (exact mass), `logp, hbd, hba, tpsa, heavy_atom_count, error_message`.

### `diagnose_smiles(smiles, try_repair=True) -> SmilesDiagnosis`

Located, categorised failure diagnosis with optional mechanical repair. See [SMILES diagnosis](../concepts/smiles-diagnosis.md).

`SmilesDiagnosis`: `input_smiles, is_valid, canonical_smiles, problems: list[SmilesProblem], repaired_smiles, repaired_is_valid, repairs_applied: list[str]`; properties `primary_category`, `caret_line()`.

`SmilesProblem`: `category, message, position` (0-based char offset or `None`), `length, atom_index, atom_indices, suggestion`.

---

## Standardisation

### `standardize_smiles(smiles, steps=None) -> StandardizeResult`

`steps` is a subset of `STANDARDIZE_STEPS` (order is fixed regardless of listing) or `None`/`['all']`. Raises `ValueError` on an unknown step.

`StandardizeResult`: `input_smiles, output_smiles, steps_applied, step_results: list[StepResult], changed, error`; property `success`. `StepResult`: `step, input_smiles, output_smiles, changed, note`.

### `enumerate_tautomers(smiles, max_tautomers=1000) -> TautomerResult`

`TautomerResult`: `input_smiles, tautomers: list[str], num_tautomers, canonical_tautomer, error`.

### `analyze_stereochemistry(smiles) -> StereoResult`

`StereoResult`: `input_smiles, chiral_centers: list[{atom_idx, config}], has_unassigned, error`.

### `get_iupac_name(smiles, use_online=False, cache_db_path=None) -> IUPACResult`

`IUPACResult`: `input_smiles, canonical_smiles, molecular_formula, molecular_weight, inchi, inchikey, iupac_name, iupac_name_source` (`pubchem` | `cache` | `None`), `error`.

---

## Molecular identity and library comparison

### `compute_identity(smiles) -> IdentityKeys`

Keys at every level of `IDENTITY_LEVELS`. `parent` is the largest fragment after neutralisation; `tautomer`, `nostereo` and `skeleton` are RDKit `RegistrationHash` layers of the parent. See [Molecular identity](../concepts/molecular-identity.md).

`IdentityKeys`: `input_smiles, is_valid, exact, parent, tautomer, nostereo, skeleton, formula, n_fragments, had_charge, has_stereo, error`; method `key(level)`.

### `group_by_identity(smiles, level='parent', keys=None) -> IdentityReport`

Pass precomputed `keys` (same order) to avoid recomputation across levels.

`IdentityReport`: `level, n_records, n_valid, n_groups` (distinct compounds at this level), `n_collapsed` (records sharing a key with an earlier record), `n_groups_by_level: dict`, `groups: list[IdentityGroup]` (multi-member only, largest first), `keys: list[IdentityKeys]`, `error`.

`IdentityGroup`: `level, key, size, indices, smiles, distinct_exact, differs_by: list[str]` — any of `salt, counter-ion or charge form` · `tautomer` · `stereochemistry` · `constitution`.

### `strictest_shared_level(a, b) -> str | None` · `describe_difference(a, b) -> str`

For two `IdentityKeys`: the strictest level at which they agree, and a plain-language description — `identical`, `salt, counter-ion or charge form`, `tautomer`, `stereochemistry`, `stereochemistry and tautomer`, `constitution (same formula only)`, `different compounds`, optionally suffixed `; also salt, counter-ion or charge form`.

### `diff_libraries(smiles_a, smiles_b, level='parent', keys_a=None, keys_b=None, include_unchanged=False) -> LibraryDiff`

`LibraryDiff`: `level, n_a, n_b, n_valid_a, n_valid_b, n_keys_a, n_keys_b, n_added, n_removed, n_unchanged, n_changed, changes_by_kind: dict[str,int], multiplicity_changes, entries: list[DiffEntry], error`; property `jaccard`.

`DiffEntry`: `status` (`added` | `removed` | `unchanged` | `changed`), `key, smiles_a, smiles_b, count_a, count_b, change`. Entries are ordered removed → added → changed → unchanged.

---

## SMARTS auditing

### `audit_smarts(patterns, library=None, library_source='user-supplied', checks=None, breadth_threshold=0.10, preparations=None, dead_sample=2500) -> SmartsAuditResult`

`patterns` is a sequence of SMARTS strings or `(smarts, name, rule_set)` tuples (as from `load_patterns`). `library` is a sequence of RDKit `Mol`; `None` loads the bundled set. `checks` is a subset of `AUDIT_CHECKS` (`compile` always runs). `preparations[0]` is the default preparation used for breadth, dead and redundancy. See [SMARTS auditing](../concepts/smarts-auditing.md).

`SmartsAuditResult`: `n_patterns, n_molecules, library_source, checks_run, breadth_threshold, patterns: list[PatternAudit], sensitivity: SensitivitySummary | None, elapsed_seconds, error`; count properties `n_unparseable, n_needs_explicit_h, n_over_broad, n_dead, n_dead_never_matching_atom, n_duplicates, n_equivalent, n_subsumed, n_clean`; method `to_rows()` for CSV.

`PatternAudit`: identity (`index, smarts, name, rule_set`); compile (`parses, parse_error, n_query_atoms, requires_explicit_h, has_recursive_smarts`); breadth (`n_hits, hit_fraction, over_broad`); dead (`never_fires, dead_verdict` ∈ {`rare combination`, `never-matching atom`, `fires only with <prep>`}, `never_matching_atoms`); redundancy (`duplicate_of, equivalent_to, subsumed_by` as indices into `patterns`); sensitivity (`hits_by_preparation, preparation_sensitive`); properties `flags: list[str]`, `clean`.

`SensitivitySummary`: `preparations, n_molecules, compounds_flagged, total_hits, patterns_firing` (each `dict[prep,int]`), `verdict_flips: dict[prep,int]` (vs the default preparation), `n_sensitive_patterns`.

### `explain_smarts(smarts, library=None, preparations=None, n_examples=5, sample=2500) -> SmartsExplanation`

`SmartsExplanation`: `smarts, parses, parse_error, normalized_smarts, n_query_atoms, n_query_bonds, requires_explicit_h, has_recursive_smarts, atoms: list[AtomExplanation], hits_by_preparation, n_molecules, example_matches, never_matching_atoms`; property `verdict`.

`AtomExplanation`: `atom_index, query` (the single-atom SMARTS), `n_matching_molecules, is_hydrogen`.

### `load_patterns(path) -> list[tuple[str, str|None, str|None]]`

CSV/TSV/XLSX with a `smarts` (or `pattern`/`query`/`smarts_pattern`) column and optional `description`/`name`/`rule_id`/`id`/`alert`/`label` and `rule_set_name`/`rule_set`/`set`/`catalog`/`source`; or text with one SMARTS per line (optional trailing name, `#` comments). Raises `FileNotFoundError`, `ValueError`.

### `load_reference_library(path=None, max_molecules=None) -> (list[Mol], str)`

`None` loads the bundled ChEMBL-derived set (see [Data and licensing](../project/data-and-licensing.md)); otherwise any supported compound file.

### `prepare_molecule(mol, preparation) -> Mol`

Copy of `mol` as `implicit-h` (unchanged), `explicit-h` (`AddHs`) or `kekule` (Kekulé bonds, aromatic flags cleared). Raises `ValueError` on an unknown preparation.

---

### `subsumes(a, b) -> ProofResult` · `equivalent(a, b) -> ProofResult` · `satisfiable(smarts) -> SatisfiabilityResult`

Static proofs. `ProofResult`: `relation, a, b, status ('proven' | 'no witness' | 'budget exceeded' | 'not analysable: …'), proven, witness (B atom → A atom), witness_reverse, status_reverse, universe_size`. `SatisfiabilityResult`: `smarts, status, satisfiable, atom_index, atom_expression, universe_size`.

### `prove_catalogue(patterns, max_container_atoms=40, step_budget=200000, progress_callback=None) -> CatalogueProof`

`patterns` is a list of `(smarts, name)`. `CatalogueProof`: `n_patterns, universe_size, universe_domains, n_analysable, n_not_analysable, not_analysable_reasons, n_unsatisfiable, n_proven_redundant, n_proven_equivalent, n_no_witness, n_pairs_tested, n_pairs_budget_exceeded, patterns: list[PatternProof]`; `PatternProof`: `index, name, smarts, status, reason, proven_subsumed_by, proven_equivalent_to, witness, n_undecided_pairs, n_budget_exceeded`. `audit_smarts(..., checks=[..., "proof"])` fills `proven_subsumed_by`, `proven_equivalent_to`, `proof_status`, `proof_reason` on each `PatternAudit`.

### `diff_smarts(source_a, source_b, library=None, preparation='implicit-h', max_molecules=None, n_examples=3) -> SmartsDiffResult`

Semantic diff of two catalogues. Sources are pattern files or `rdkit:<NAME>` (`RDKIT_CATALOGS`). `SmartsDiffResult`: `source_a, source_b, library_source, n_molecules, preparation, n_patterns_a, n_patterns_b, n_paired, paired_by, text_counts, semantic_counts, flagged_a, flagged_b, flagged_both, flagged_only_a, flagged_only_b, verdict_changes, verdict_change_fraction, patterns: list[PatternDiff]`. `PatternDiff`: `key, paired_by, a, b (PatternSide: name, smarts, rule_set, parses, n_hits), text_status, semantic_status, hits_a, hits_b, gained, lost, jaccard, examples_gained, examples_lost`.

---

## Screening and search

### `apply_filters(smiles, rules=None, preparation='implicit-h') -> FilterResult`

`rules` ⊆ `{lipinski, veber, ghose, egan, ro3, pains, qed}` or `None`/`['all']`. `preparation` affects only PAINS matching and is recorded.

`FilterResult`: `smiles, molecular_weight, logp, hbd, hba, tpsa, rotatable_bonds, heavy_atom_count, molar_refractivity, qed_score`, per-rule `RuleResult` fields `lipinski, veber, ghose, egan, ro3, pains`, `passes_all, preparation, error`. `RuleResult`: `passed, details`.

### `substructure_search(query, library, is_smarts=True, preparation='implicit-h') -> list[SubstructureHit]`

`SubstructureHit`: `smiles, matched, match_indices, preparation`.

### `compute_similarity(query_smiles, library, fp_type='ecfp4', n_bits=2048, threshold=0.0, top_n=None) -> list[SimilarityResult]`

`SimilarityResult`: `rank, query, hit, similarity, fingerprint_type`.

### `compute_fingerprint(smiles, fp_type='ecfp4', n_bits=2048) -> FingerprintResult | list[FingerprintResult]`

`fp_type='all'` returns a list, one per type. `FingerprintResult`: `smiles, fingerprint_type, n_bits, n_on_bits, density, bit_string, hex_string`.

### `extract_scaffold(smiles) -> ScaffoldResult`

`ScaffoldResult`: `input_smiles, scaffold_smiles` (empty for acyclic), `note, error`.

### `rgroup_decomposition(core_smarts, smiles_list) -> list[RGroupResult]`

`RGroupResult`: `input_smiles, decomposition: dict[str,str]` (`Core`, `R1`, `R2` …), `is_matched, error`.

---

## Structures

### `generate_structure(smiles, output_path, format='sdf', dimension='3d', force=False, title=None) -> str`

Writes a file; returns a status string (`Generated`, `Skipped (File already exists)`, or an error description). Formats `sdf`, `mol`, `pdb`.

### `generate_conformers(smiles, num_conformers=50, output_sdf=None) -> ConformerResult`

`ConformerResult`: `input_smiles, num_generated, error`.

### `validate_reaction(smiles) -> ReactionResult` · `map_atoms(smiles) -> AtomMapResult` · `augment_smiles(smiles, num_augmentations=5) -> AugmentResult`

`ReactionResult`: `input_smiles, is_valid, num_reactants, num_products, num_agents, error`. `AtomMapResult`: `input_smiles, mapped_smiles, error`. `AugmentResult`: `input_smiles, augmented_smiles: list[str], error`.

---

## Databases (network)

### `resolve(query, sources=None, query_type='auto', unichem=True, timeout=40.0, use_cache=True) -> ResolveResult`

Query several databases in parallel under one time budget and reconcile the answers. `sources` is a list drawn from `PROVIDERS` (`pubchem`, `chembl`, `chebi`, `kegg`) or `['all']`; default `DEFAULT_SOURCES = ['pubchem', 'chembl', 'chebi']`. `query_type` is one of `auto`, `name`, `smiles`, `inchikey`, `id`. Raises `ValueError` on an unknown source or type; network failures never raise — they appear as `error` outcomes.

`ResolveResult`: `query, query_type, sources, records: list[CompoundRecord], outcomes: list[SourceOutcome], merged: CompoundRecord | None, cross_refs: dict[str, str], consensus_inchikey, agreement ('agree' | 'disagree' | 'unknown'), found`; `by_source(name)` returns that source's record or `None`.

`ResolveResult` also carries `agreement_level` (strictest identity level shared by every returned structure: `exact`, `parent`, `tautomer`, `nostereo`, `skeleton`, `formula`, `none` for different compounds, or `None` when fewer than two structures are available), `disagreement` (plain-language summary), `pairwise: list[PairDifference]` and `corrections: list[Correction]` (`source, text_hit_id, text_hit_name, text_hit_inchikey, difference, corrected, corrected_id`).

`SourceOutcome`: `source, status ('found' | 'not found' | 'error'), seconds, error`.

### `concordance(queries, sources=None, query_type='auto', timeout=40.0, use_cache=True, progress_callback=None) -> ConcordanceReport`

Resolve every query (UniChem off) and tabulate agreement. `ConcordanceReport`: `sources, n_queries, n_found_any, n_found_all, agreement_counts, level_counts, disagreement_classes, found_by_source, errors_by_source, corrections_by_source, correction_classes_by_source, rows: list[ConcordanceRow]`; `ConcordanceRow`: `query, n_found, agreement, agreement_level, disagreement, consensus_inchikey, ids, names, inchikeys, corrections, errors`.

`CompoundRecord` (one schema for every database): `source, source_id, name, synonyms, smiles, inchi, inchikey, formula, molecular_weight, monoisotopic_mass, charge, xlogp, tpsa, hbd, hba, rotatable_bonds, cross_refs, url, extra`.

### `get_provider(name) -> Provider`

Cached provider instance. `Provider.lookup(query, query_type='auto', deadline=None) -> CompoundRecord | None`; `by_id`, `by_name`, `by_smiles`, `by_inchikey`; `looks_like_id(query)`. Raises `ProviderError` on HTTP failure. Implementations: `PubChemProvider`, `ChEMBLProvider`, `ChEBIProvider`, `KEGGProvider` in `chemlitmus.providers`.

### `unichem_xrefs(inchikey, timeout=15.0) -> dict[str, str]`

Identifiers for one InChIKey across the [UniChem](https://www.ebi.ac.uk/unichem/) sources, keyed by a short source name (`chembl`, `drugbank`, `chebi`, `pubchem`, `cas`, `hmdb`, `kegg`, `pdbe`, `surechembl`, …; unknown sources as `unichem_src_<n>`).

### PubChem-specific

### `lookup(query, search_type='auto', use_cache=True) -> PubChemCompound | None`

### `lookup_by_name(name, use_cache=True) -> PubChemCompound | None`

### `lookup_file(input_file, output_file=None, output_format='csv', remove_duplicates=True, progress_callback=None) -> list[PubChemCompound]`

### `download_structure(cid, format='sdf', dimension='3d', output_dir='structures', force=False) -> str`

`PubChemCompound`: `input_query, cid, title, iupac_name, molecular_formula, molecular_weight, canonical_smiles, isomeric_smiles, inchikey, inchi, xlogp, hbond_donor_count, hbond_acceptor_count`.

Lower-level: `PubChemClient` (rate-limited HTTP client) and `DatabaseManager` (SQLite cache) are exported for advanced use; see their docstrings.

---

## Utilities

### `parse_compounds_file(path) -> list[str]`

The input reader every batch command uses. `from chemlitmus.utils.parsers import parse_compounds_file`.

### `export_results(results, output_path, fmt)`

Write `PubChemCompound` lists to `csv` | `xlsx` | `json`. `from chemlitmus.utils.export import export_results`.

### `chemlitmus.config.settings`

The resolved `Settings` object; see [Configuration](../getting-started/configuration.md).
