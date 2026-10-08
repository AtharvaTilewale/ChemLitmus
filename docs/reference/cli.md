# CLI reference

```
chemlitmus [--version] COMMAND [OPTIONS] [ARGS]
```

Every command has `--help`. Conventions used throughout:

- **Single vs batch.** Most commands take a SMILES as a positional argument *or* a file via `--file`/`-f` (`-i` on the older `filter`, `fingerprint`, `similar` and `download`). Batch mode prints a summary table and writes full results with `--output`/`-o`.
- **Input files.** `.csv`, `.tsv`, `.xlsx`, `.smi`, `.sdf`, `.txt`. A column named `smiles`, `canonical_smiles`, `structure` or `compound` (any case) is used; otherwise the first column. In `.smi`, the first whitespace-delimited token per line is the SMILES.
- **Exit codes.** `0` success · `1` usage or runtime error · `2` input judged invalid (`validate`, `diagnose` in single mode) or a `--chiral-flag` failure. Batch commands exit `0` even when some records fail; per-record status is in the output.
- **Quoting.** Single-quote any SMARTS containing `$`, and any SMILES containing `(`, `[` or `#` in shells that treat them specially.

---

## Validation and diagnosis

### `validate`

Validate SMILES offline and report canonical form plus basic properties.

```
chemlitmus validate [SMILES] [--file PATH] [--output PATH] [--quiet]
```

| Option | Description |
|---|---|
| `SMILES` | single SMILES to validate |
| `--file`, `-f` | batch input |
| `--output`, `-o` | save batch results to CSV |
| `--quiet`, `-q` | single mode: print only the canonical SMILES (nothing on failure) |

Single mode exits `2` on an invalid SMILES. Whitespace inside a SMILES is treated as invalid (RDKit would silently parse only the first token). Batch CSV columns: `input_smiles, is_valid, canonical_smiles, molecular_formula, molecular_weight, logp, hbd, hba, tpsa, heavy_atom_count, error_message`.

### `diagnose`

Explain why a SMILES fails to parse, with character positions and suggested fixes; attempt safe mechanical repairs.

```
chemlitmus diagnose [SMILES] [--file PATH] [--output PATH] [--no-repair] [--only-invalid | --all]
```

| Option | Description |
|---|---|
| `--file`, `-f` | batch input |
| `--output`, `-o` | per-record diagnoses to CSV |
| `--no-repair` | do not attempt repairs |
| `--only-invalid` / `--all` | batch: report only invalid records (default) or every record |

Checks, in order: `characters`, `brackets`, `parentheses`, `rings`, `syntax`, `valence`, `aromaticity`. Single mode exits `2` on an invalid SMILES. CSV columns: `input_smiles, is_valid, canonical_smiles, primary_category, n_problems, problems, positions, suggestions, repaired_smiles, repaired_is_valid, repairs_applied`. See [SMILES diagnosis](../concepts/smiles-diagnosis.md).

---

## Standardisation and identity

### `standardize`

Clean a SMILES through salt stripping, neutralisation, tautomer canonicalisation and canonical output.

```
chemlitmus standardize [SMILES] [--steps LIST] [--file PATH] [--output PATH] [--show-diff]
```

| Option | Description |
|---|---|
| `--steps`, `-s` | comma-separated subset of `fragment, neutralize, tautomer, canonical`, or `all` (default). Steps always run in that order regardless of how they are listed |
| `--show-diff` | single mode: show the SMILES after each step |
| `--file`, `--output` | batch CSV: `input_smiles, output_smiles, changed, error` |

### `identity`

Compute layered identity keys for one molecule, or group a collection by identity level.

```
chemlitmus identity [SMILES] [--file PATH] [--level LEVEL] [--output PATH] [--show N]
```

| Option | Description |
|---|---|
| `--level`, `-l` | `exact` · `parent` (default) · `tautomer` · `nostereo` · `skeleton` · `formula` |
| `--output`, `-o` | per-record keys at every level plus `group_id_<level>`, `group_size_<level>`, `group_varies_by` |
| `--show` | multi-member groups to list in the terminal (default 15) |

See [Molecular identity](../concepts/molecular-identity.md).

### `diff`

Structure-aware comparison of two compound collections.

```
chemlitmus diff FILE_A FILE_B [--level LEVEL] [--output PATH] [--json PATH] [--include-unchanged] [--show N]
```

| Option | Description |
|---|---|
| `FILE_A`, `FILE_B` | reference (older) and comparison (newer) collections |
| `--level`, `-l` | identity level for matching (default `parent`) |
| `--output`, `-o` | CSV: `status, level, key, smiles_a, smiles_b, count_a, count_b, change` |
| `--json` | complete result including summary counts |
| `--include-unchanged` | also list unchanged compounds |
| `--show` | entries per category in the terminal (default 10) |

### `tautomers`

```
chemlitmus tautomers [SMILES] [--max N] [--file PATH] [--output PATH]
```

Enumerate tautomers (RDKit `TautomerEnumerator`); the canonical one is marked. `--max`/`-m` caps enumeration (default 1000). Batch output has one row per tautomer: `input_smiles, tautomer_smiles, is_canonical, error`. Exits `1` on an invalid single SMILES.

### `stereo`

```
chemlitmus stereo SMILES [--chiral-flag]
```

List stereocentres with CIP labels (`?` = unassigned). `--chiral-flag` exits `1` if any centre is unassigned — a CI gate.

### `iupacname`

```
chemlitmus iupacname [SMILES] [--online] [--file PATH] [--output PATH]
```

InChI, InChIKey, formula and exact mass offline; `--online` adds the PubChem preferred IUPAC name (cached after first fetch).

---

## Pattern quality control

### `smartsaudit`

Audit a SMARTS pattern set against a reference molecule population, or explain one pattern.

```
chemlitmus smartsaudit [PATTERNS] [--explain SMARTS] [--library PATH] [--checks LIST]
                       [--breadth-threshold F] [--max-molecules N] [--output PATH] [--json PATH] [--show N]
```

| Option | Description |
|---|---|
| `PATTERNS` | CSV/TSV/XLSX with a `smarts` column (optional `description`/`name`, `rule_set_name`/`rule_set`), or text with one SMARTS per line |
| `--explain`, `-e` | decompose a single SMARTS instead of auditing a file |
| `--library`, `-l` | reference molecules; default is the bundled ChEMBL-derived set (9,272) |
| `--checks`, `-c` | subset of `compile, breadth, dead, redundancy, sensitivity` (default all; `compile` always runs) |
| `--breadth-threshold` | hit fraction above which a pattern is over-broad (default 0.10) |
| `--max-molecules` | use only the first N reference molecules (faster, coarser) |
| `--output`, `-o` | per-pattern table CSV |
| `--json` | full result including the sensitivity summary |
| `--show` | flagged patterns to list (default 15) |

Per-pattern CSV columns: `index, name, rule_set, smarts, parses, parse_error, n_query_atoms, requires_explicit_h, has_recursive_smarts, n_hits, hit_fraction, over_broad, never_fires, dead_verdict, never_matching_atoms, duplicate_of, equivalent_to, subsumed_by, preparation_sensitive, flags, hits_implicit-h, hits_explicit-h, hits_kekule`. See [Auditing alert sets](../guides/audit-alert-sets.md) and [SMARTS auditing](../concepts/smarts-auditing.md).

---

## Screening and search

### `filter`

Drug-likeness and ADMET rules with PAINS alerts.

```
chemlitmus filter [SMILES] [--file PATH] [--rules LIST] [--fail] [--qed-min F] [--output PATH] [--prep PREP]
```

| Option | Description |
|---|---|
| `--file`, `-i` | batch input |
| `--rules`, `-r` | subset of `lipinski, veber, ghose, egan, ro3, pains, qed`, or `all` (default). `qed` is reported but never fails a compound |
| `--fail` | keep only compounds that **fail** (e.g. to extract PAINS hits) |
| `--qed-min` | batch: drop compounds below this QED |
| `--prep` | molecule preparation for substructure-based rules: `implicit-h` (default) · `explicit-h` · `kekule`. Recorded in output |
| `--output`, `-o` | CSV with properties, per-rule pass/fail, `passes_all`, `preparation` |

Rule definitions: Lipinski MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10 · Veber RotB ≤ 10, TPSA ≤ 140 · Ghose 160 ≤ MW ≤ 480, −0.4 ≤ LogP ≤ 5.6, 20 ≤ heavy atoms ≤ 70, 40 ≤ MR ≤ 130 · Egan TPSA ≤ 131.6, LogP ≤ 5.88 · Ro3 MW ≤ 300, LogP ≤ 3, HBD ≤ 3, HBA ≤ 3 · PAINS: RDKit `FilterCatalog` PAINS A/B/C.

### `substructure`

```
chemlitmus substructure QUERY --file PATH [--output PATH] [--smiles-query] [--prep PREP]
```

SMARTS query by default; `--smiles-query` treats it as an exact SMILES fragment. Output: `smiles, match_indices, preparation`.

### `similar`

```
chemlitmus similar QUERY --file PATH [--threshold F] [--top N] [--fp-type T] [--bits N] [--output PATH]
```

Tanimoto ranking. Defaults: threshold 0.5, top 10, `ecfp4`, 2048 bits. Output: `rank, query, hit, similarity, fingerprint_type`.

### `fingerprint`

```
chemlitmus fingerprint [SMILES] [--file PATH] [--type T] [--bits N] [--output PATH]
```

`--type`: `ecfp4, ecfp6, fcfp4, maccs, rdkit, atompair, torsion, all`. MACCS is fixed at 167 bits. Batch output is hex-encoded, one column per type.

### `scaffold`

```
chemlitmus scaffold [SMILES] [--file PATH] [--output PATH]
```

Murcko scaffold. Acyclic molecules are reported as `acyclic` (empty scaffold). CSV: `smiles, scaffold, error`.

### `rgroup`

```
chemlitmus rgroup [CORE] [--core SMARTS] [--smiles LIST] [--file PATH]
```

R-group decomposition against a core with labelled attachment points (`[*:1]`, `[*:2]` …). The core may be given positionally or with `--core`/`-c`. Molecules without the core are listed with a status, not dropped.

---

## Structures and 3D

### `download`

Fetch structures from PubChem, or generate them locally from SMILES.

```
chemlitmus download [QUERY] [--file PATH] [--format F] [--3d | --2d] [--output-dir DIR] [--force] [--gen MODE]
```

| Option | Description |
|---|---|
| `QUERY` | CID, SMILES or name |
| `--file`, `-i` | batch |
| `--format`, `-f` | `sdf` (default), `mol`, `pdb`, `png` (PubChem only) |
| `--3d` / `--2d` | coordinate dimension (default 3D) |
| `--output-dir`, `-o` | default `structures` |
| `--force` | overwrite existing files (default skips them) |
| `--gen`, `-g` | `all` — generate every structure locally, **no network**; `missing` — PubChem first, local generation when absent |

Local generation: RDKit ETKDG v3 (random-coordinate fallback) + MMFF94 (UFF fallback); 2D via `Compute2DCoords`.

### `conformers`

```
chemlitmus conformers SMILES [--num N] [--output PATH]
```

ETKDG v3 + MMFF94 conformer ensemble to a multi-model SDF. `--num`/`--num-conformers`/`-n` default 50; `--output` default `conformers.sdf`.

### `reaction`

```
chemlitmus reaction "REACTANTS>AGENTS>PRODUCTS"
```

Parse a reaction SMILES/SMIRKS; report reactant, agent and product template counts. Exits `1` if unparseable.

### `atommap`

```
chemlitmus atommap SMILES
```

Assign sequential atom-map numbers to every atom.

### `augment`

```
chemlitmus augment SMILES [--num N]
```

Randomised, non-canonical SMILES of the same molecule (default 5) for machine-learning data augmentation.

---

## PubChem

### `lookup`

```
chemlitmus lookup QUERY [--type T] [--no-cache]
```

`--type`/`-t`: `auto` (default), `cid`, `smiles`, `name`, `inchikey`. Unknown values are rejected before any network call. See [PubChem](../guides/pubchem.md).

### `batch`

```
chemlitmus batch INPUT_FILE [--output PATH] [--format csv|xlsx|json] [--keep-duplicates]
```

Multithreaded, rate-limited batch lookup with a `<output>_report.log`.

---

## Utilities

### `status`

Show version, configuration source, cache/data/log directories and database state.

### `init`

Create the cache, data and log directories and the SQLite database.

### `update`

```
chemlitmus update [--check] [--yes]
```

Check for a newer release on PyPI or GitHub and optionally install it.
