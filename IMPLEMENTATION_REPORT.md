# Implementation report

Work carried out against `CHEMLITMUS_IMPLEMENTATION_BRIEF.md` (brief written against commit
`a6ee9af`, package version 1.0.0). This document records, per feature: purpose, what changed,
affected interfaces, scientific assumptions, example usage, verification, and remaining limits.

Released as **1.1.0**. 347 offline tests pass; `ruff check chemlitmus tests --select F,E9` is
clean; `mkdocs build --strict` succeeds; the package builds.

---

## 1. Verification of the brief's P0 findings

Every source finding in §4 was checked against the working tree before any change. All five were
**confirmed**; one needed a narrower statement than the brief gave.

| Finding | Verdict | Evidence at `a6ee9af` |
|---|---|---|
| §4.1 `fcfp4` is ECFP4 | **Confirmed** | `cheminfo.py` returned `GetMorganGenerator(radius=2, fpSize=n_bits)` for both `ecfp4` and `fcfp4`; byte-identical fingerprints. |
| §4.2 inconsistent structure validation | **Confirmed** | `compute_identity` rejected internal whitespace; `standardize_smiles`, `enumerate_tautomers`, `analyze_stereochemistry`, `extract_scaffold`, `get_iupac_name`, `compute_fingerprint`, `compute_similarity` and `substructure_search` all called `Chem.MolFromSmiles` directly, so `"CC O"` became ethane in eight entry points. |
| §4.3 failures reported as chemical negatives | **Confirmed** | `_match_vector` returned an all-false vector on any matching exception; `prepare_molecule` returned the *unprepared* molecule when Kekulisation failed, so a Kekulé audit silently mixed preparations. |
| §4.4 parser discards records | **Confirmed** | `parse_compounds_file` returned a list of strings, dropped every other column, removed blanks, skipped unparseable SDF entries, and guessed the first column. |
| §4.5 repairs presented without caution | **Confirmed, partially mitigated already** | The model exposed `repaired_smiles` / `repaired_is_valid` / `repairs_applied` and the docs warned that repairs are mechanical, but there was no status field, no policy gate, and the CLI printed "Repaired:". |
| §4.6 scientific descriptions | **Confirmed for one claim** | `docs/concepts/molecular-identity.md` already stated the branching relationship correctly ("`tautomer` and `nostereo` are siblings"), but the worked example immediately below it read consecutive counts as exclusive causes ("578 are salt variants; a further 187 are tautomer variants"). That paragraph was wrong and is rewritten. |

One claim in the repository that the brief did not name was found to be wrong during §12 work and
is corrected: see [§8](#8-corrected-quantitative-claim).

---

## 2. P0 — corrected semantics

### 2.1 Fingerprint semantics

**Purpose.** A fingerprint name must describe the computation performed.

**Changes.** `fcfp4` now builds a Morgan generator with
`atomInvariantsGenerator=GetMorganFeatureAtomInvGen()` — pharmacophoric (donor / acceptor /
aromatic / halogen / basic / acidic) atom invariants, which is what FCFP means.
`FingerprintResult` gained `algorithm`, `radius`, `atom_invariants`, `representation` and
`chirality`; `FP_DESCRIPTIONS` documents every supported name.

**Interfaces.** `compute_fingerprint`, `compute_similarity` (both CLI and API) — additive fields;
**the bits produced for `fcfp4` change**, so any stored `fcfp4` fingerprints must be recomputed.

**Verification.** Four molecules compared bit-for-bit against RDKit's own feature-Morgan and
connectivity-Morgan generators. The aggregate test asserts ECFP and FCFP differ for most of a
panel, not for every individual molecule (they legitimately coincide for simple structures).

### 2.2 One strict parser

**Purpose.** RDKit's `MolFromSmiles("CC O")` returns ethane with "O" as a title. Any structure
field with internal whitespace must be an error, not a different molecule.

**Changes.** `chemlitmus.core.smiles.mol_from_smiles()` is now the single entry point; it raises
`SmilesParseError` for empty or whitespace-bearing input, except for a CXSMILES extension block
(`CCO |$;;R1$|`), which is parsed with `allowCXSMILES`. `split_smiles_field()` handles `.smi`
lines, which *define* a name after whitespace. `standardize`, `tautomers`, `stereo`, `scaffold`,
`iupacname` and every `cheminfo` entry point route through it.

**Assumption.** A structure *column* in a table holds one SMILES. A `.smi` *line* may hold a name.
CXSMILES is the only sanctioned internal whitespace.

**Verification.** Parametrised tests assert `"CC O"` produces an error in every consumer and that
CXSMILES still parses.

### 2.3 Failures are not chemical negatives

**Purpose.** An incomplete audit must not claim that a rule never fires or that a compound passes.

**Changes.** `PreparedLibrary` builds the substructure library *and* records which molecules
reached the requested preparation; a molecule that cannot be Kekulised is **unevaluated**, not a
non-hit, and is excluded from hit masks and verdict-flip counts. `prepare_molecule_status()`
returns `(mol, "ok")` or `(None, reason)`; `prepare_molecule(..., strict=True)` raises
`PreparationError`. Matching errors become `PatternAudit.match_error` (hit counts stay `None`),
the `match-failed` flag, and `SmartsAuditResult.n_match_failures`; such patterns are never called
dead. `SmartsAuditResult.preparation_status` reports evaluated / failed per preparation, and
`hit_fraction` now divides by *evaluated* molecules. `smartsdiff` excludes unprepared molecules
and reports `n_unevaluated`. The dataset audit exposes `processing_complete`.

**Verification.** Tests monkey-patch a Kekulisation failure and a matcher exception and assert the
accounting, the absence of a dead verdict, and the `match-failed` flag.

### 2.4 Record-preserving ingestion

**Purpose.** Every downstream count must reconcile with the input.

**Changes.** New `chemlitmus.core.records`: `Record` (stable `record_id`, `position`, `source_id`,
original `structure`, all other `fields`, `status`, `issues`) and `RecordSet` (counts per status,
detected `roles`, `notes`, `reconcile()`). CSV/TSV/TXT/XLSX/SMI/SDF adapters. SDF properties are
preserved and unreadable records kept with `SDF_PARSE_FAILED`. Ambiguous column mapping raises
`SchemaError` rather than guessing; `ambiguous="first"` opts into the old behaviour.
`parse_compounds_file()` remains as a documented compatibility wrapper.

**Assumption.** Source ids need not be unique; internal ids are position-based and always unique.

**Verification.** Tests cover duplicate source ids, empty fields, whitespace fields, invalid
SMILES, `.smi` comments and names, SDF properties plus a corrupted record, XLSX, explicit roles,
and both ambiguity paths. `n_total == sum(statuses) == len(records)` is asserted on every path.

### 2.5 Repairs are candidates

**Changes.** `SmilesDiagnosis.repair_status` ∈ {`not attempted`, `not needed`, `candidate`,
`failed`}. A candidate parses; it is not claimed to be the intended molecule and carries no
confidence number. `ChemicalPolicy.repair.mode` (`none` / `report` / `apply`) gates whether
candidates are produced or written; the original record is always retained. The CLI prints
"Candidate repair (parses; not verified to be the intended molecule)".

### 2.6 Corrected descriptions

The identity report paragraph now states that the rows answer independent questions, that the
counts are not additive because `tautomer` and `nostereo` are siblings below `parent`, that
differences between consecutive rows are not causes, and that equal formulas identify a formula
group rather than a compound. `formula` is excluded from `ChemicalPolicy.identity_level` by
validation.

---

## 3. P0 — policy, provenance and reproducibility

**`ChemicalPolicy`** (validated, hashed, serialisable) covers representations, fragment selection,
neutralisation, tautomer / stereo / isotope handling, identity level, match preparation and
chirality, alert sets, fingerprint and similarity threshold, repair mode and per-code severity
overrides. Presets: `conservative` (keep everything, identity at `exact`, no repairs) and
`parent`. Multiple substantial organic components are **flagged for review**
(`FRAG_MULTIPLE_ORGANIC`), not silently resolved by size.

**`Transformation`** records `before`, `after`, operation, parameters, `changed_text`,
`changed_identity` and — when identity changed — the relation in identity-level language. Atom
correspondence is **not** claimed; the brief permits reporting it only where supported, and no
supported case is implemented.

**`RunManifest`** records tool / Python / RDKit versions, platform, input and output checksums
(file names only — no private paths), policy hash and name, rule-catalogue identity, settings and
output schema version.

**Limit.** Reproducibility holds *within a supported environment*. RDKit updates can change
aromaticity perception, standardisation and descriptors; the manifest makes such a change
detectable but does not prevent it.

---

## 4. P1 — unified dataset audit

**Purpose.** One command, one offline pass, every record accounted for.

```bash
chemlitmus audit dataset.csv --output-dir audit_results --config policy.json \
    --split-column split --endpoint-column pIC50 --units-column units --clean --fail-on-severity error
```

**Changes.** `chemlitmus.core.dataset_audit` with 28 codes in `ISSUE_CATALOGUE`, each mapped to a
default severity, a meaning and a suggested action. Per record: parse status, standardisation with
provenance, component / charge / isotope / stereo / element flags, descriptors, Murcko scaffold,
alerts with matched atom indices, identity keys at all six levels, group membership, and duplicate
codes per level (reported separately, strictest level per pair only, so they are not double
counted). Dataset level: identity groups with `differs_by` and the splits they span, descriptor
summaries, alert prevalence with its denominator, and dataset-level warnings kept distinct from
record-level issues.

`chemlitmus.core.report` writes `records.csv` (all original columns prefixed `in:`), `issues.csv`,
`identity_groups.csv`, `summary.json`, `audit.json`, `policy.json`, `manifest.json`, an offline
`report.html` and, with `--clean`, `clean.csv` plus `exclusions.csv` listing every exclusion and
its reason. The HTML escapes all user text, renders depictions robustly (invalid structures become
a label), highlights matched atoms, shows denominators beside every metric, and provides a
client-side filter over the issue table.

**Gates.** `GatePolicy` → exit **0** pass, **1** configuration error, **3** policy violation,
**4** partial processing. No universal opaque score is produced.

**Verification.** A 13-record fixture exercising salts, duplicates, whitespace, empty fields,
invalid SMILES, censored measurements, stereo variants, isotopes and a two-component mixture; the
test asserts the exact issue-code set per record, exact accounting, `clean + exclusions == input`,
manifest contents, HTML escaping, and all six gate outcomes.

**Limit.** Descriptor rules and alert matches describe structure. They do not establish toxicity,
ADMET behaviour, activity or synthesisability, and the report says so.

---

## 5. P1 — leakage and label conflicts

**Leakage** (`chemlitmus.core.leakage`) reports exact / parent / tautomer / nostereo / skeleton
overlap as separate classes with their own counts, example pairs and fractions. **The classes nest
and are not additive**, which the report, the CLI and the HTML all state. Scaffold overlap and
nearest-neighbour Tanimoto (fingerprint and threshold named) are reported as *relatedness*;
acyclic molecules are excluded from scaffold overlap with their count given. Formula matches are
reported and explicitly not called identity leakage. Within-split duplicates are counted per
split. A temporal check runs when dates exist and says so when they do not.

**Conflicts** (`chemlitmus.core.labels`) compares records only within the same identity group *and*
the same endpoint context (target / source / assay columns). Molar units convert to nM with the
conversion string retained and a pX value derived; the spread is measured in log10 units
(tolerance 1.0 = ten-fold) or in absolute units when a single non-molar unit is shared. Censored
values (`<`, `>`, `<=`, `>=`, `~`, ranges) are kept as bounds and never averaged; missing units
are flagged and excluded from comparison. The default action is review.

**Verification.** Curated salt / tautomer / stereoisomer cases land in the correct classes;
constitutional isomers sharing a formula do **not** become identity leakage; unit conversion,
censoring, context separation and classification conflicts are all asserted.

---

## 6. P2 — splits and generation

**`make_splits`** supports `random`, `identity`, `scaffold`, `temporal` and `source`. Groups are
indivisible; **each acyclic molecule forms its own scaffold group** rather than collapsing into one
giant group. Seeds reproduce assignments. The report gives achieved fractions, the group-size
distribution including the largest group as a share of the data, endpoint balance per split, and
every unsatisfiable constraint as an explicit message. The result is verified with the same
leakage implementation the audit uses.

**`evaluate_generated`** reports validity over *all attempts*, uniqueness over *valid outputs* at
each identity level, and novelty over *unique valid outputs* against the named reference sets only.
Repaired candidates are evaluated as a separate population. When no valid output exists, metrics
are listed under `undefined` rather than fabricated as zero. No learned-model distribution metrics
are included (they would require model identity and environment provenance the brief rightly
demands; none is implemented, and none is claimed).

**Verification.** Salt / stereo novelty is asserted to vary with the identity policy
(1.0 at `exact`, 0.5 at `parent`, 0.0 at `nostereo` for the same inputs); impossible split
constraints are asserted to be reported rather than resolved by breaking a group; the empty-valid
case is asserted to return `undefined`.

---

## 7. P2 — SMARTS evidence, provenance and holdout

`smartsaudit` now records:

* **evidence labels** — `exact text duplicate`, `identical observed hit set`,
  `observed hit-set containment`, `static proven containment`, `static proven equivalence`,
  `not observed in this reference`, `unsupported`, `undecided` — so an observed claim about one
  panel is never read as a universal one;
* **match semantics** — preparation, hydrogen handling, chirality (ignored, as RDKit does by
  default), aromaticity model and RDKit version;
* **catalogue metadata** — name / version / source / licence supplied by the user, plus a content
  hash and the pattern count **per rule set**;
* **reference-panel composition** — size, median heavy atoms, ring and charge fractions, element
  census, and a content hash;
* **example matches with matched atoms**, and examples of molecules that match only under a
  non-default preparation;
* **`--holdout`** libraries — additional panels whose hits name the patterns that never fired on
  the main panel, so "never fires" is never extrapolated into impossibility.

The static prover's existing states (`no witness`, `budget exceeded`, `not analysable`) remain
undecided and are mapped onto the evidence labels unchanged. Local atom/bond satisfiability is
still documented as *not* a proof that a whole molecule is realisable.

**Not implemented** (and not claimed): precision/recall against curated positive/negative labels
(no independently curated label set is available to this work — manufacturing one would be the
very kind of unsupported claim the brief warns against), cross-version RDKit regression runs, and
expansion of the prover's supported primitives.

---

## 8. Corrected quantitative claim

The README stated: *"the identical PAINS catalogue flags 77.6%, 90.5% or 88.5% of compounds"*.
Re-measured with the new versioned recipe
(`benchmarks/alert_preparation_sensitivity.py`, catalogue SHA-256 `92f29232…`, bundled
9,272-molecule panel):

| Rule set | Patterns | implicit-H / explicit-H / Kekulé | Flips vs implicit-H |
|---|---|---|---|
| PAINS | 481 | 3.8% / 5.0% / 2.1% | 1.3% / 4.2% |
| SureChEMBL | 166 | 26.2% / 68.7% / 41.2% | 42.5% / 24.6% |
| Inpharmatica | 91 | 26.2% / 26.8% / 48.1% | 0.5% / 28.2% |
| all eight sets | 1,251 | 77.6% / 90.5% / 88.5% | 16.5% / 16.6% |

The 77.6/90.5/88.5 figures are for the **union of eight published sets**, not PAINS. README and
`docs/concepts/molecule-preparation.md` now give the per-set table and state the denominator; the
underlying scientific point (preparation changes verdicts, and nobody records it) is *stronger*
per set, not weaker: SureChEMBL flips 42.5% of the library.

---

## 9. Worked examples

### Research data curator — reconcile a messy table

```bash
chemlitmus audit curation.csv --endpoint-column value --units-column units \
    --relation-column relation --target-column assay --clean -o audit_results
```

On an 8-record fixture: 6 ok, 1 empty, 1 invalid; 4 distinct compounds at `parent` with 2 records
collapsing; errors — 1 `PARSE_WHITESPACE` (with the first token shown as evidence); warnings —
`LABEL_CONFLICT` ×2, `STD_CHANGED_IDENTITY`, `DUP_PARENT`, `PARSE_EMPTY`, `DUP_EXACT`; info —
`FRAG_MULTIPLE`, `ELEMENT_UNUSUAL`, `REPAIR_CANDIDATE`, `RELATION_CENSORED`, `STEREO_UNASSIGNED`.
`clean.csv` and `exclusions.csv` sum back to 8.

### ML engineer — split and leakage audit

```bash
chemlitmus split dataset.csv --strategy scaffold -f train=0.8,test=0.2 --seed 0 -o split.csv
chemlitmus audit dataset.csv --split-column split --no-split-overlap --fail-on-severity error
```

The split report names every group it could not divide; the audit's gate fails the run (exit 3) if
any evaluation record matches a training record at the policy identity level.

### Screening researcher — alert review

```bash
chemlitmus smartsaudit alerts.csv --checks all,proof --holdout vendor_deck.smi \
    --catalogue-name "in-house alerts" --catalogue-version 2026.1 -o audit.csv
```

Each pattern carries its evidence label, example matches with matched atoms, and — for patterns
that fire only on the holdout — proof that a "never fires" verdict was panel-specific.

### Generative chemist — evaluate samples

```bash
chemlitmus generated samples.smi -r train.smi --constraints 'mw=0:500,logp=-1:5' --repair -o per_molecule.csv
```

Validity over all attempts, uniqueness over valid outputs, novelty over unique valid outputs
against `train.smi` only, with repaired candidates reported as their own population.

---

## 9a. Feature matrix against the brief

Every numbered requirement of the brief, with its state. "Partial" and "not implemented" are
stated as such; nothing below is marked complete on the strength of a placeholder.

| Brief | Requirement | State | Where |
|---|---|---|---|
| §4.1 | FCFP4 uses feature invariants; fingerprint provenance recorded | **Implemented** | `cheminfo.FP_DESCRIPTIONS`, `FingerprintResult` |
| §4.2 | One shared structure validator; `.smi` names vs structure columns; defined CXSMILES behaviour | **Implemented** | `smiles.mol_from_smiles`, `split_smiles_field`, `is_cxsmiles` |
| §4.3 | Execution/preparation status; failures are not negatives | **Implemented** | `smartsaudit.PreparedLibrary`, `PreparationStatus`, `match_error`, `processing_complete` |
| §4.4 | Record-preserving ingestion; explicit column mapping; compatibility wrapper | **Implemented** | `core/records.py`, `utils/parsers.parse_compounds_file` |
| §4.5 | Repairs as candidates with edits, rationale, status; policy gate; no fabricated confidence | **Implemented** | `SmilesDiagnosis.repair_status`, `RepairPolicy` |
| §4.6 | Corrected identity / formula / policy / proof descriptions | **Implemented** | `docs/concepts/molecular-identity.md`, `smarts-proofs.md` |
| §5 | Typed record model, stable ids, preserved failures, SDF properties | **Implemented** | `Record`, `RecordSet` |
| §5 | CSV/TSV/SMI/SDF/XLSX adapters | **Implemented** | `core/records.py` |
| §5 | Parquet as an optional extension | **Implemented** (extra `chemlitmus[parquet]`) | `records._read_parquet` |
| §5 | Versioned, validated, hashed policy with presets | **Implemented** | `core/policy.py` |
| §5 | Mixtures flagged rather than resolved by size; scope status for unsupported chemistry | **Implemented** | `FRAG_MULTIPLE_ORGANIC`, `ELEMENT_UNUSUAL`, `status='unsupported'` |
| §5 | Transformation provenance; atom correspondence where supported | **Implemented / deliberately absent** | `Transformation`; atom correspondence is **not** claimed anywhere |
| §5 | Run manifest with checksums, versions, policy hash, seeds | **Implemented** | `core/manifest.py` |
| §6 | One audit command and Python entry point; required checks | **Implemented** | `audit_dataset`, `chemlitmus audit` |
| §6 | Stable codes, severity, record ids, evidence, action; dataset vs molecule findings | **Implemented** | `ISSUE_CATALOGUE` (28 codes), `AuditSummary.dataset_warnings` |
| §6 | Annotated records, quarantine table, groups, summary JSON, manifest, offline HTML; clean export with stated policy | **Implemented** | `core/report.py` |
| §6 | No universal opaque score | **Implemented** (none exists) | — |
| §7 | Five overlap classes as distinct evidence; within-split duplicates | **Implemented** | `core/leakage.py` |
| §7 | Scaffold overlap and configurable nearest-neighbour similarity, parameters reported | **Implemented** | `PairReport.neighbour_fingerprint/threshold` |
| §7 | Formula matches never called identity leakage | **Implemented** | `PairReport.formula_matches` |
| §7 | Optional temporal/source constraints; unavailable metadata identified | **Implemented** | `temporal_violations`, `temporal_note` |
| §7 | Measured contamination distinguished from unknown pretraining contamination | **Implemented** | `LeakageReport.scope` |
| §8 | Conflicts within identity groups; endpoint context required | **Implemented** | `core/labels.py`, `context_fields` |
| §8 | Unit conversion with provenance; explicit log conversion; originals preserved | **Implemented** | `Measurement.conversion`, `value_nm`, `log_value` |
| §8 | Censored values treated as bounds, never averaged | **Implemented** | `Measurement.censored`, excluded from spread |
| §8 | Technical replicates separated from biological variability when metadata allows | **Implemented** | `ReplicateSplit` (`source_field`) |
| §8 | Default to review, not strongest activity | **Implemented** | `ConflictGroup.suggested_action = "review"` |
| §9 | Offline HTML with plots, searchable tables, depictions, highlighted atoms, escaped text | **Implemented** | `report.render_html` |
| §9 | JSON authoritative, CSV per record, versioned schemas | **Implemented** | `audit.json`, `OUTPUT_SCHEMA_VERSION` |
| §9 | Severity thresholds and explicit gates; stable exit behaviour | **Implemented** | `GatePolicy`, exits 0/1/3/4 |
| §9 | GitHub Actions and Python gate examples | **Implemented** | `docs/guides/audit-a-dataset.md` |
| §10 | Random / identity / scaffold / temporal / source splits; groups indivisible; seeds | **Implemented** | `core/splits.py` |
| §10 | Achieved fractions, balance, group sizes, post-split leakage; conflicts exposed | **Implemented** | `SplitReport` |
| §10 | Acyclic molecules not pooled into one scaffold group | **Implemented** | `_scaffold_key` |
| §11 | Validity / uniqueness / novelty with declared denominators; repaired evaluated separately | **Implemented** | `core/generation.py` |
| §11 | Scaffold diversity, neighbours, descriptors, alerts, constraints | **Implemented** | `GenerationReport` |
| §11 | Undefined metrics on empty valid output | **Implemented** | `GenerationReport.undefined` |
| §11 | Uncertainty intervals only with a justified sampling model | **Deliberately absent** — none is justified here, so none is reported | — |
| §11 | Learned-model distribution metrics, optional, with model provenance | **Not implemented** | would need model identity + environment capture |
| §12 | Catalogue metadata (source, version, licence, reference, domain) | **Implemented** | `CatalogueMetadata` |
| §12 | Example molecules, atom highlights, preparation-flip and version-flip examples | **Implemented** | `example_matches`, `example_match_atoms`, `preparation_flip_examples`; version flips via `smartsdiff` |
| §12 | Shared match-semantics configuration recorded in results | **Implemented** | `MatchSemantics` |
| §12 | Precise evidence labels | **Implemented** | `EVIDENCE_LABELS` (8) |
| §12 | Holdout evaluation; no extrapolation of zero hits | **Implemented** | `--holdout`, `HoldoutResult` |
| §12 | Reference-panel composition summaries; configurable panels | **Implemented** | `ReferencePanel`, `--library` |
| §12 | Reviewable cleanup proposals preserving published provenance | **Implemented** | `core/cleanup.py`, `smartsaudit --cleanup` |
| §12 | Undecided prover states preserved; satisfiability is local | **Implemented** | unchanged from 1.0.0, documented |
| §12 | Independently checkable witnesses, adversarial cases | **Implemented** | witnesses returned; 20 directional cases in `tests/test_smartsproof.py` |
| §12 | Precision/recall against curated labels | **Not implemented** | no independently curated label set available |
| §12 | Cross-version RDKit regressions | **Not implemented** | needs a multi-version matrix |
| §12 | Versioned benchmark recipe; no union-as-PAINS claims | **Implemented** | `benchmarks/alert_preparation_sensitivity.py`; claim corrected |
| §13 | Streaming/chunked ingestion, resumable jobs, on-disk identity index | **Not implemented** | see §10 |
| §13 | Packed/sparse SMARTS matrices | **Not implemented** | dense boolean matrix retained |
| §13 | Throughput and peak-memory benchmarks | **Not implemented** | so no capacity is claimed anywhere |
| §13 | Provider retrieval timestamps, cache staleness, field-level provenance | **Implemented** | `CompoundRecord.retrieved_at/from_cache/cached_at/field_provenance` |
| §13 | Source record versions where available | **Partial** | field exists (`source_version`); no provider populates it — none of the four exposes one |
| §13 | Cache snapshot export/import | **Not implemented** | |
| §13 | Core audits run without network | **Implemented** | audit/leakage/conflicts/split/generated make no network call |
| §13 | Optional REST/MCP adapter, report viewer | **Not implemented** | deferred as the brief directs |
| §14 | Curated regression corpus with reasons and provenance | **Implemented** | `tests/test_brief_coverage.py::CORPUS` (18 cases) |
| §14 | Verification against direct RDKit reference operations | **Implemented** | FCFP4 vs RDKit generators; identity vs `RegistrationHash` |
| §14 | Standardisation idempotence, identity invariance, group-preserving splits, accounting, CLI/API parity | **Implemented** | `tests/test_brief_coverage.py` |
| §14 | Deterministic offline suite in CI; larger benchmarks separate; mocked providers | **Implemented** | 376 offline tests; `benchmarks/`; provider tests mock HTTP |
| §14 | Documented supported RDKit range | **Implemented** | `pyproject.toml` (`RDKit>=2023.09`), stated in `docs/getting-started/installation.md` |
| §14 | Published benchmark inputs, checksums, configurations, outputs, limitations | **Implemented** | benchmark script prints and stores the catalogue checksum |
| §14 | Controlled "cleaning improves model performance" comparison | **Not attempted** | and no such claim is made |
| §15 | Shared typed analysis layer; thin adapters; additive fields; optional dependency groups | **Implemented** | `core/*` + `cli/commands/*`; extras `dev`, `docs`, `parquet` |
| §16 | `IMPLEMENTATION_REPORT.md` with purpose, changes, interfaces, assumptions, usage, verification, limits | **Implemented** | this document |

## 10. Remaining work (explicitly unfinished)

| Item | Why not done |
|---|---|
| Bounded-memory / streaming execution (§13) | The audit holds the record set and a dense pattern matrix in memory. Chunked ingestion, an on-disk identity index and resumable jobs are designed for but not implemented; no throughput or peak-memory benchmark is published, so no capacity is claimed. |
| REST / MCP adapter and report viewer (§13) | Deferred until the core workflows settle, as the brief directs. |
| Provider **source record versions** (§13) | The `source_version` field exists but no provider populates it: none of PubChem, ChEMBL, ChEBI or KEGG exposes a per-record version in its API responses. Retrieval timestamps, cache staleness and field-level provenance *are* implemented. |
| Cache snapshot export/import (§13) | Licensing varies per source; not attempted. |
| SMARTS precision/recall and cross-version RDKit regressions (§12) | Requires an independently curated label set and a multi-version test matrix; neither exists here, and fabricating either would be an unsupported claim. |
| Controlled "cleaning improves model performance" comparison (§14) | Not attempted. No such claim is made anywhere in the code or documentation. |

## 11. Compatibility notes

* **`fcfp4` bits change.** Any stored `fcfp4` fingerprints or similarity results must be recomputed.
* **Structure fields with internal whitespace now error** where they previously parsed as the first
  token. This is the intended fix; datasets that relied on the old behaviour will surface records
  that were silently wrong.
* **`parse_compounds_file`** keeps its signature and return type.
* All new result fields are additive; no existing field was removed or repurposed.
* New optional dependencies (`openpyxl`, `matplotlib`) are declared in `pyproject.toml`.
