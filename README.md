<p align="center">
  <img src="https://raw.githubusercontent.com/AtharvaTilewale/ChemLitmus/main/docs/assets/chemlitmus-logo-white.png" alt="ChemLitmus" width="560" />
</p>

<h3 align="center">Quality control for small-molecule data</h3>

<p align="center">
  Validate, standardise, deduplicate, audit and compare chemical datasets, structural-alert catalogues and database records — before anything is modelled, screened or published.
</p>

<p align="center">
  <a href="https://github.com/AtharvaTilewale/ChemLitmus/actions/workflows/ci.yaml"><img src="https://github.com/AtharvaTilewale/ChemLitmus/actions/workflows/ci.yaml/badge.svg" alt="CI"></a>
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue" alt="Python 3.10–3.13"></a>
  <a href="https://www.rdkit.org/"><img src="https://img.shields.io/badge/RDKit-%E2%89%A5%202023.09-green" alt="RDKit >= 2023.09"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-yellow.svg" alt="MIT"></a>
  <a href="https://chemlitmus.readthedocs.io"><img src="https://img.shields.io/badge/docs-chemlitmus.readthedocs.io-informational" alt="Documentation"></a>
</p>

<p align="center">
  <a href="#installation">Installation</a> ·
  <a href="#quick-start">Quick start</a> ·
  <a href="#what-chemlitmus-does">Features</a> ·
  <a href="#typical-workflows">Workflows</a> ·
  <a href="#python-api">Python API</a> ·
  <a href="#commands-at-a-glance">Commands</a> ·
  <a href="#documentation">Documentation</a> ·
  <a href="#citing">Citing</a>
</p>

---

## Overview

Chemical datasets are rarely as clean as they look. SMILES that RDKit silently truncates at a space, the same compound entered as three salt forms, a test set that shares parents with the training set, a pIC50 of 11.2 sitting next to neighbours at 8, an alert catalogue whose verdicts change depending on how hydrogens are written — none of these are caught by `MolFromSmiles` returning a molecule.

**ChemLitmus** is a command-line tool and Python library that finds these problems and reports them as *evidence* a scientist can act on: every record is accounted for, every finding carries a stable code with the affected records and a suggested action, every number states its denominator, and nothing is silently dropped, averaged or "fixed".

- **35 commands**, all offline on RDKit except the explicitly network-backed database lookups
- **Typed results** (Pydantic models) for every operation; CSV / JSON / HTML outputs for every batch command
- **A declared chemical policy** — fragment handling, charges, tautomers, stereo, identity level, molecule preparation, alert sets — that is validated, hashed and written into every output
- **Stable exit codes** so the same checks run unchanged in CI

## Installation

```bash
pip install chemlitmus
```

| Alternative | Command |
|---|---|
| Isolated CLI only | `pipx install chemlitmus` |
| uv | `uv tool install chemlitmus` (CLI) · `uv add chemlitmus` (dependency) |
| Parquet input support | `pip install "chemlitmus[parquet]"` |
| From source, for development | `git clone https://github.com/AtharvaTilewale/ChemLitmus.git && cd ChemLitmus && pip install -e ".[dev]"` |

Requirements: Python 3.10 or newer; RDKit 2023.09 or newer is installed automatically. Linux, macOS and Windows are tested in CI. Every dependency ships as a binary wheel — no compiler needed.

```bash
chemlitmus --version      # v1.0.0
chemlitmus --help
```

## Quick start

Every file below is in the repository's [`examples/`](examples/) directory, with a README describing what each one deliberately gets wrong.

**Why is this SMILES broken, and what would fix it?**

```bash
$ chemlitmus diagnose "c1cncc1"
  c1cncc1
  ^
  [aromaticity] Can't kekulize mol.  Unkekulized atoms: 0 1 2 3 4  (position 0)
      → A pyrrole-type nitrogen must carry its hydrogen: write it as [nH]
  Candidate repair (parses; not verified to be the intended molecule): c1cc[nH]c1
```

**How many distinct compounds are in this file, and what varies between the duplicates?**

```bash
$ chemlitmus identity --file examples/demo.csv
┃ Level    ┃ Distinct compounds ┃ Collapsed records ┃
│ exact    │                  3 │                 0 │
│ parent < │                  2 │                 1 │
...
│    2 │ salt, counter-ion or charge form │ CC(=O)Oc1ccccc1C(=O)O | CC(=O)Oc1ccccc1C(=O)O.[Na+] │
```

**Audit a dataset in one pass — parsing, standardisation, duplicates, alerts, leakage, label conflicts:**

```bash
$ chemlitmus audit examples/demo.csv --endpoint-column pIC50 --split-column split -o audit/
│ Records                     │     6 │ ok 3 · empty 1 · invalid 2 · unsupported 0 · error 0 │
│ Distinct compounds (parent) │     2 │ 1 records collapse onto an earlier one               │
│ Issues: error               │     3 │ SPLIT_OVERLAP ×1, PARSE_WHITESPACE ×1, PARSE_INVALID ×1 │
│ Leakage test→train          │     1 │ exact 0 · parent 1 · … (nested, not additive)        │
│ Label conflicts             │     1 │ 1 groups compared (quantitative)                     │
Outputs: audit/  (records.csv, issues.csv, identity_groups.csv, summary.json, audit.json, policy.json, report.html, manifest.json)
```

**Is this structural-alert catalogue sound?**

```bash
$ chemlitmus smartsaudit examples/alerts.csv --checks all,proof
# unparseable, dead, over-broad, duplicate, equivalent and subsumed patterns;
# verdict sensitivity to molecule preparation; static containment proofs
```

**Does one query agree across databases?**

```bash
$ chemlitmus resolve "levothyroxine" --sources pubchem,chembl,chebi,kegg
# per-source structures, the kind of any disagreement (salt form, tautomer, stereo), merged record
```

## What ChemLitmus does

| Area | What you get | Commands |
|---|---|---|
| **Dataset audit** | One offline pass: parse status, standardisation provenance, identity groups, structural flags, descriptors, alerts, split leakage, label conflicts. 32 stable issue codes with severity, evidence and action. HTML + JSON + CSV outputs, run manifest, optional clean export listing every exclusion, CI gates. | `audit` |
| **Leakage and conflicts** | Train/test overlap as *nested* evidence classes (exact → parent → tautomer / nostereo → skeleton), scaffold and similarity relatedness reported separately; contradictory labels among records of the same compound with molar units converted and censored values kept as bounds. | `leakage` `conflicts` |
| **Activity cliffs** | Matched molecular pairs and near neighbours whose labels are *proven* to differ; label outliers whose every neighbour disagrees with them — SAR or error, reported for review. | `cliffs` |
| **ML data workflows** | Reproducible splits that never divide an identity or scaffold group; generated-molecule evaluation with an explicit denominator for every metric. | `split` `generated` |
| **Validation and diagnosis** | Located, explained SMILES failures with candidate repairs; catches records RDKit would silently truncate at whitespace. | `validate` `diagnose` |
| **Standardisation and identity** | Salt stripping, neutralisation, tautomer canonicalisation; six nested identity levels on RDKit `RegistrationHash`; tautomer enumeration, stereo analysis, InChI / IUPAC identifiers. | `standardize` `identity` `tautomers` `stereo` `iupacname` |
| **Collection comparison** | Structure-aware diff of two files: added, removed, unchanged, and *changed with the reason* (salt form, tautomer, stereo). | `diff` |
| **SMARTS catalogue QC** | Defect audit of alert sets; static containment *proofs* with checkable witnesses (no reference library needed); behavioural diff between two catalogue versions, or a file against RDKit's built-in PAINS / Brenk / NIH. | `smartsaudit` `smartsproof` `smartsdiff` |
| **Screening and search** | Lipinski / Veber / Ghose / Egan / Ro3, PAINS and QED with the molecule preparation declared and recorded; substructure, similarity, fingerprints, scaffolds, R-group decomposition. | `filter` `substructure` `similar` `fingerprint` `scaffold` `rgroup` |
| **Structures and 3D** | 2D / 3D generation (ETKDG + MMFF94), conformer ensembles, reaction SMILES, atom mapping, SMILES augmentation. | `download --gen all` `conformers` `reaction` `atommap` `augment` |
| **Databases** | One query across PubChem, ChEMBL, ChEBI and KEGG with UniChem cross-references and a merged record; name → structure concordance statistics for a whole list; cached lookups and batch retrieval. | `resolve` `concordance` `lookup` `batch` `download` |

## Typical workflows

<details>
<summary><b>Audit a dataset before training a model</b></summary>

```bash
# 1. One pass over everything, with CI gates
chemlitmus audit examples/bioactivity.csv --endpoint-column standard_value --units-column standard_units \
    --split-column split -o audit/ --no-split-overlap --max-invalid-fraction 0.01

# 2. Look at leakage between the splits on its own
chemlitmus leakage examples/bioactivity.csv --split-column split --json leakage.json

# 3. Contradictory measurements of the same compound (censored values stay bounds)
chemlitmus conflicts examples/bioactivity.csv -e standard_value --units-column standard_units \
    --relation-column standard_relation --context assay_id --source-column document_id

# 4. Near-identical compounds whose labels disagree; suspect labels for review
chemlitmus cliffs examples/bioactivity.csv -e standard_value --units-column standard_units \
    --relation-column standard_relation -o pairs.csv --outliers outliers.csv

# 5. Rebuild the split so no identity or scaffold group is divided
chemlitmus split examples/bioactivity.csv -s scaffold -f train=0.8,test=0.2 --seed 0 -o split.csv
```

Exit codes from `audit`: `0` pass · `1` configuration error · `3` a gate was violated · `4` processing incomplete.
Guide: [Audit a dataset](docs/guides/audit-a-dataset.md) · [Splits and generated molecules](docs/guides/splits-and-generation.md)
</details>

<details>
<summary><b>Clean and deduplicate a compound library</b></summary>

```bash
chemlitmus diagnose --file examples/library.csv -o diagnosed.csv   # every failure located and explained
chemlitmus standardize --file examples/library.csv -o std.csv       # salts, charges, tautomers — with provenance
chemlitmus identity --file std.csv --level parent -o groups.csv # what is actually the same compound
chemlitmus audit std.csv -o audit/ --clean                      # clean export + the list of every exclusion
```

Guide: [Clean a compound library](docs/guides/clean-a-library.md)
</details>

<details>
<summary><b>Audit and compare structural-alert catalogues</b></summary>

```bash
chemlitmus smartsaudit examples/alerts.csv --checks all,proof -o audit.csv    # defects + sensitivity + proofs
chemlitmus smartsproof examples/alerts.csv                                     # containment proofs, no library needed
chemlitmus smartsdiff examples/alerts.csv examples/alerts_v2.csv               # what changed, in molecules not text
chemlitmus smartsdiff examples/alerts.csv rdkit:PAINS                          # a file against RDKit's built-in set
chemlitmus filter --file examples/library.csv --rules pains --prep explicit-h  # screen with the preparation recorded
```

Guide: [Audit a structural-alert set](docs/guides/audit-alert-sets.md)
</details>

<details>
<summary><b>Look compounds up across databases</b></summary>

```bash
chemlitmus resolve "aspirin"                                  # PubChem + ChEMBL + ChEBI + KEGG, reconciled
chemlitmus resolve CHEMBL25 --json aspirin.json
chemlitmus concordance --file examples/names.txt -o concordance.csv   # how often do sources agree on a structure?
chemlitmus lookup 2244 --type cid
chemlitmus batch examples/smiles.txt -o metadata.csv
```

Guide: [Database lookups](docs/guides/databases.md)
</details>

## Python API

Everything the CLI does is available as functions returning typed models.

```python
from chemlitmus import (
    ChemicalPolicy, read_records, audit_dataset, activity_cliffs,
    diagnose_smiles, group_by_identity, diff_libraries, audit_smarts, apply_filters,
)

# Diagnose, don't just reject
d = diagnose_smiles("CC O")
d.is_valid, d.repaired_smiles, d.repair_status
# (False, 'CCO', 'candidate')   # a candidate parses; it is not claimed to be the intended molecule

# Audit a table under a declared policy
policy = ChemicalPolicy.preset("parent")
records = read_records("examples/bioactivity.csv", roles={"endpoint": "pchembl_value", "split": "split"})
audit = audit_dataset(records, policy)
audit.summary.issues_by_code          # {'SPLIT_OVERLAP': 1, 'DUP_PARENT': 1, 'PARSE_WHITESPACE': 1, ...}
audit.summary.processing_complete     # True — every record reached a terminal status

# Activity cliffs and label outliers
recs = [dict(r.fields, record_id=r.record_id, smiles=r.parsed_smiles) for r in records.ok_records()]
cliffs = activity_cliffs(recs, policy, endpoint_field="pchembl_value")
cliffs.n_cliffs, cliffs.n_undetermined, [o.record_ids for o in cliffs.outliers]

# Identity, comparison, catalogue QC, screening
group_by_identity(smiles_list, level="parent").n_groups
diff_libraries(old_smiles, new_smiles).changes_by_kind      # {'salt, counter-ion or charge form': 3, ...}
audit_smarts(smarts_patterns).n_subsumed                      # measured on the bundled 9,272-molecule panel
apply_filters("CC(=O)Oc1ccccc1C(=O)O", ["pains"], preparation="explicit-h").preparation   # 'explicit-h'
```

Full reference: [Python API](docs/reference/python-api.md).

## Design principles

These are the rules every command follows; they are what distinguish a finding from a guess.

| Principle | In practice |
|---|---|
| **Nothing is dropped silently** | Every input record has a status (`ok`, `empty`, `invalid`, `unsupported`, `error`) and the counts reconcile exactly. Ambiguous column mapping is an error, never a guess. |
| **Evidence, not verdicts** | Findings are issue codes with affected records, evidence and a suggested action. No opaque quality score. Default action on a conflict is *review*, not "keep the strongest". |
| **Every number names its denominator** | Validity over all attempts, uniqueness over valid outputs, novelty over unique valid outputs; alert rates per rule set, never a union figure attributed to one set. |
| **Censored values are bounds** | `> 10 µM` is never averaged with `5 µM`. A pair of measurements is a conflict or a cliff only when the bounds *prove* the difference; otherwise it is undetermined. |
| **A failed computation is not a chemical negative** | A pattern that raised is "unknown", not "zero hits"; a molecule that could not be prepared is "unevaluated", not "clean". |
| **Repairs are candidates** | A candidate repair parses. It carries no confidence number and is applied only under an explicit repair policy, with the original retained. |
| **Proven is kept apart from observed** | A SMARTS containment proof holds for all molecules; a hit count holds for the panel it was measured on. The reports say which is which. |
| **The chemistry is declared** | Fragment, charge, tautomer, stereo, identity level, molecule preparation, alert sets and fingerprint are a versioned, hashed policy recorded in every output. |

## Why molecule preparation is recorded

On 9,272 ChEMBL molecules, the same alert catalogue gives different answers depending only on how the molecules were prepared — implicit hydrogens, explicit hydrogens or Kekulé bonds — a choice most tools neither expose nor record. Measured with [`benchmarks/alert_preparation_sensitivity.py`](benchmarks/alert_preparation_sensitivity.py) on the 1,251 ChEMBL structural alerts (`rd_filters` redistribution, eight published sets):

| Rule set | Patterns | Compounds flagged (implicit-H / explicit-H / Kekulé) | Verdict flips vs implicit-H |
|---|---|---|---|
| SureChEMBL | 166 | 26.2% / 68.7% / 41.2% | 42.5% / 24.6% |
| Inpharmatica | 91 | 26.2% / 26.8% / 48.1% | 0.5% / 28.2% |
| Glaxo | 55 | 12.7% / 10.6% / 36.9% | 2.1% / 24.7% |
| PAINS | 481 | 3.8% / 5.0% / 2.1% | 1.3% / 4.2% |
| **All eight sets** | 1,251 | 77.6% / 90.5% / 88.5% | 16.5% / 16.6% |

42.5% of compounds change their SureChEMBL verdict between two defensible preparations, and only 19 of the 481 PAINS patterns fire at all on this library. ChemLitmus therefore reports alert figures per set with their denominator, takes the preparation as an explicit `--prep` option, and writes it into every output row.

## Commands at a glance

| Group | Commands |
|---|---|
| Dataset quality | `audit` · `leakage` · `conflicts` · `cliffs` · `split` · `generated` |
| Validation | `validate` · `diagnose` |
| Standardisation and identity | `standardize` · `identity` · `tautomers` · `stereo` · `iupacname` · `diff` |
| Pattern quality control | `smartsaudit` · `smartsproof` · `smartsdiff` |
| Screening and search | `filter` · `substructure` · `similar` · `fingerprint` · `scaffold` · `rgroup` |
| Structures and 3D | `download` · `conformers` · `reaction` · `atommap` · `augment` |
| Databases | `resolve` · `concordance` · `lookup` · `batch` |
| Utilities | `init` · `status` · `update` |

`chemlitmus <command> --help` documents every option. Conventions shared by all commands:

- **Input** — SMILES on the command line, or `--file` / a positional dataset in CSV, TSV, XLSX, SMI, SDF or Parquet. Column roles (structure, id, endpoint, split, units, …) are detected from headers or given explicitly.
- **Output** — a readable terminal report; `-o` for CSV; `--json` for the full typed result. Batch commands never stop at the first bad record.
- **Exit codes** — `0` success · `1` usage or runtime error · `2` input judged invalid (single-molecule `validate` / `diagnose`) · `3` an `audit` gate violated · `4` `audit` processing incomplete.
- **Configuration** — `--config policy.json` (or `-c`; see [`examples/policy.json`](examples/policy.json)) selects the chemical policy; `chemlitmus status` shows paths, cache and the resolved defaults.

## Documentation

Full documentation: **[chemlitmus.readthedocs.io](https://chemlitmus.readthedocs.io)** (built from [`docs/`](docs/)).

| | |
|---|---|
| Getting started | [Installation](docs/getting-started/installation.md) · [Quickstart](docs/getting-started/quickstart.md) · [Example data](examples/README.md) · [Configuration](docs/getting-started/configuration.md) |
| Guides | [Audit a dataset](docs/guides/audit-a-dataset.md) · [Splits and generated molecules](docs/guides/splits-and-generation.md) · [Clean a library](docs/guides/clean-a-library.md) · [Audit an alert set](docs/guides/audit-alert-sets.md) · [Compare collections](docs/guides/compare-libraries.md) · [Search and structures](docs/guides/search-and-structures.md) · [Database lookups](docs/guides/databases.md) |
| Reference | [CLI](docs/reference/cli.md) · [Python API](docs/reference/python-api.md) |
| Concepts | [Records, policy and provenance](docs/concepts/records-and-policy.md) · [Molecular identity](docs/concepts/molecular-identity.md) · [Activity cliffs](docs/concepts/activity-cliffs.md) · [Molecule preparation](docs/concepts/molecule-preparation.md) · [SMARTS auditing](docs/concepts/smarts-auditing.md) · [SMARTS containment proofs](docs/concepts/smarts-proofs.md) · [SMILES diagnosis](docs/concepts/smiles-diagnosis.md) |
| Project | [Changelog](CHANGELOG.md) · [Contributing](docs/project/contributing.md) · [Citing](docs/project/citing.md) · [Data and licensing](docs/project/data-and-licensing.md) |

## Project status

Version 1.0.0. Tested on Python 3.10–3.13 across Linux, macOS and Windows with every push; 393 offline tests plus live-database integration tests; documentation built in strict mode. Deliberately out of scope for this release: streaming execution for datasets larger than memory, activity prediction, and any claim that cleaning improves a particular model — the tool reports what the data contains and leaves the modelling conclusions to you.

## Contributing

Bug reports, feature requests and pull requests are welcome through [GitHub issues](https://github.com/AtharvaTilewale/ChemLitmus/issues). The development setup, test commands and documentation build are described in [Contributing](docs/project/contributing.md). In short:

```bash
pip install -e ".[dev]"
ruff check chemlitmus tests --select F,E9
pytest -q -m "not integration"
mkdocs build --strict
```

## Citing

If ChemLitmus contributes to published work, please cite it — see [`CITATION.cff`](CITATION.cff) and [Citing](docs/project/citing.md) for BibTeX. Results that use the bundled reference library should also cite ChEMBL (Zdrazil *et al.*, *Nucleic Acids Res.* 2024), and results that use the ChEMBL structural alerts should cite their original publications as listed in [Data and licensing](docs/project/data-and-licensing.md).

## Licence

Source code is released under the [MIT Licence](LICENSE).

The bundled reference library (`chemlitmus/data/reference_library.smi.gz`, 9,272 molecules) is derived from ChEMBL 37 and redistributed under [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/); it is **not** covered by the MIT licence. Details in [Data and licensing](docs/project/data-and-licensing.md).

## Acknowledgements

ChemLitmus is built on [RDKit](https://www.rdkit.org/). The structural-alert benchmark uses the ChEMBL alert collection as redistributed by Pat Walters' [`rd_filters`](https://github.com/PatWalters/rd_filters). Database commands use the public web services of [PubChem](https://pubchem.ncbi.nlm.nih.gov/), [ChEMBL](https://www.ebi.ac.uk/chembl/), [ChEBI](https://www.ebi.ac.uk/chebi/), [KEGG](https://www.kegg.jp/) and [UniChem](https://www.ebi.ac.uk/unichem/).
