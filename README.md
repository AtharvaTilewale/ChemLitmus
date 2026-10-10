# ChemLitmus

<p align="center">
  <img src="https://raw.githubusercontent.com/AtharvaTilewale/ChemLitmus/main/docs/assets/chemlitmus-logo-white.png" alt="ChemLitmus" width="620" />
</p>

<p align="center">
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.10%2B-blue" alt="Python 3.10+"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="MIT"></a>
  <a href="https://github.com/AtharvaTilewale/ChemLitmus/actions/workflows/ci.yaml"><img src="https://github.com/AtharvaTilewale/ChemLitmus/actions/workflows/ci.yaml/badge.svg" alt="CI"></a>
</p>

**A command-line tool and Python library for checking, cleaning, comparing and auditing small-molecule data — the step before any modelling or screening can be trusted.**

```bash
pip install chemlitmus
```

```bash
chemlitmus diagnose "c1cncc1"          # where and why a SMILES is broken, with a fix
chemlitmus identity --file lib.csv     # how many distinct compounds, and what varies
chemlitmus diff old.smi new.smi        # what changed between two releases, chemically
chemlitmus smartsaudit alerts.csv      # is this PAINS/alert catalogue actually sound?
chemlitmus filter --file lib.csv --rules pains --prep explicit-h   # screen, reproducibly
```

Offline on RDKit, except the explicitly network-backed database commands. Every result is a typed Pydantic model; every batch command writes CSV; every command has a stable exit code.

## What it does

| | Commands |
|---|---|
| **Audit a dataset** — one offline pass: parsing, standardisation provenance, identity groups, structural flags, descriptors, alerts, train/test leakage, label conflicts; every record accounted for, stable issue codes with evidence and actions, HTML + JSON + CSV outputs, CI gates; activity cliffs between matched molecular pairs with censored values kept as bounds, and label outliers whose every neighbour disagrees with them | `audit` `leakage` `conflicts` `cliffs` |
| **ML data workflows** — group-aware splits that never divide an identity or scaffold group; generated-molecule evaluation with explicit denominators for validity, uniqueness and novelty | `split` `generated` |
| **Validate and diagnose** — located, explainable SMILES failures with mechanical repairs; catches records RDKit silently truncates at whitespace | `validate` `diagnose` |
| **Standardise and identify** — salt stripping, neutralisation, tautomer canonicalisation; six nested identity levels (exact → parent → tautomer/nostereo → skeleton → formula) built on RDKit `RegistrationHash` | `standardize` `identity` `tautomers` `stereo` `iupacname` |
| **Compare collections** — structure-aware diff: added, removed, unchanged, and *changed with reason* (salt form, tautomer, stereo) | `diff` |
| **Audit and compare SMARTS catalogues** — unparseable, dead, over-broad, duplicate, equivalent and subsumed patterns; sensitivity of verdicts to molecule preparation; static containment *proofs* with checkable witnesses (no library needed); semantic diff of two catalogue versions (or a file against RDKit's built-in PAINS/Brenk/NIH) with the number of molecules whose verdict changes | `smartsaudit` `smartsproof` `smartsdiff` |
| **Screen and search** — Lipinski/Veber/Ghose/Egan/Ro3, PAINS and QED with the molecule preparation declared and recorded; substructure, similarity, fingerprints, scaffolds, R-groups | `filter` `substructure` `similar` `fingerprint` `scaffold` `rgroup` |
| **Structures** — 2D/3D generation (ETKDG + MMFF94), conformer ensembles, reaction SMILES, atom maps, SMILES augmentation | `download --gen all` `conformers` `reaction` `atommap` `augment` |
| **Databases** — one query across PubChem, ChEMBL, ChEBI and KEGG; agreement check with the *kind* of disagreement (salt form, tautomer, stereo), merged record, UniChem cross-references; name→structure concordance statistics for a whole list | `resolve` `concordance` `lookup` `batch` `download` |

## Why these features

On 9,272 ChEMBL molecules, the same alert catalogue gives different answers depending only on how the molecules were prepared — implicit hydrogens, explicit hydrogens or Kekulé bonds — a choice no tool records. Measured with [`benchmarks/alert_preparation_sensitivity.py`](benchmarks/alert_preparation_sensitivity.py) on the 1,251 ChEMBL structural alerts (`rd_filters` redistribution, eight published sets):

| Rule set | Patterns | Compounds flagged (implicit-H / explicit-H / Kekulé) | Verdict flips vs implicit-H |
|---|---|---|---|
| SureChEMBL | 166 | 26.2% / 68.7% / 41.2% | 42.5% / 24.6% |
| Inpharmatica | 91 | 26.2% / 26.8% / 48.1% | 0.5% / 28.2% |
| Glaxo | 55 | 12.7% / 10.6% / 36.9% | 2.1% / 24.7% |
| PAINS | 481 | 3.8% / 5.0% / 2.1% | 1.3% / 4.2% |
| **all eight sets** | 1,251 | 77.6% / 90.5% / 88.5% | 16.5% / 16.6% |

Those are different questions with different answers: 42.5% of compounds change SureChEMBL verdict between two defensible preparations, and only 19 of the 481 PAINS patterns fire at all on this library. ChemLitmus reports the figure *per set with its denominator*, records the preparation in every output row, and separates what was proven from what was merely observed.

That is the gap ChemLitmus fills: the inputs are not clean, the filters are not right, and the only way to know is to measure. `smartsaudit` measures the catalogue; `--prep` makes the choice explicit and writes it into every output row; `diagnose` and `identity` do the same for the molecules.

## Documentation

**[chemlitmus.readthedocs.io](https://chemlitmus.readthedocs.io)**

- [Installation](docs/getting-started/installation.md) · [Quickstart](docs/getting-started/quickstart.md) · [Configuration](docs/getting-started/configuration.md)
- Guides: [Clean a library](docs/guides/clean-a-library.md) · [Audit an alert set](docs/guides/audit-alert-sets.md) · [Compare collections](docs/guides/compare-libraries.md) · [Search and structures](docs/guides/search-and-structures.md) · [Database lookups](docs/guides/databases.md)
- Reference: [CLI](docs/reference/cli.md) · [Python API](docs/reference/python-api.md)
- Concepts: [Molecular identity](docs/concepts/molecular-identity.md) · [Molecule preparation](docs/concepts/molecule-preparation.md) · [SMARTS auditing](docs/concepts/smarts-auditing.md) · [SMILES diagnosis](docs/concepts/smiles-diagnosis.md)

## Python

```python
from chemlitmus import diagnose_smiles, group_by_identity, diff_libraries, audit_smarts, apply_filters

d = diagnose_smiles("CC O")
d.is_valid, d.problems[0].message, d.repaired_smiles
# (False, "Whitespace splits the record: RDKit silently parses only 'CC' …", 'CCO')

rep = group_by_identity(smiles_list, level="parent")
rep.n_groups, rep.groups[0].differs_by
# (11902, ['salt, counter-ion or charge form'])

res = audit_smarts(patterns)                       # bundled 9,272-molecule reference set
res.n_dead, res.n_subsumed, res.sensitivity.verdict_flips

apply_filters("CC(=O)Oc1ccccc1C(=O)O", ["pains"], preparation="explicit-h").preparation
# 'explicit-h'
```

## Requirements

Python 3.10+ · RDKit ≥ 2023.09 · pandas · NumPy · Typer · Rich · Pydantic. All install as wheels; no compiler needed.

## Licence

Source code: [MIT](LICENSE). The bundled reference library (`chemlitmus/data/`) is derived from ChEMBL 37 and redistributed under [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/) — see [Data and licensing](docs/project/data-and-licensing.md). To cite, see [`CITATION.cff`](CITATION.cff).
