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
| **Validate and diagnose** — located, explainable SMILES failures with mechanical repairs; catches records RDKit silently truncates at whitespace | `validate` `diagnose` |
| **Standardise and identify** — salt stripping, neutralisation, tautomer canonicalisation; six nested identity levels (exact → parent → tautomer/nostereo → skeleton → formula) built on RDKit `RegistrationHash` | `standardize` `identity` `tautomers` `stereo` `iupacname` |
| **Compare collections** — structure-aware diff: added, removed, unchanged, and *changed with reason* (salt form, tautomer, stereo) | `diff` |
| **Audit SMARTS catalogues** — unparseable, dead, over-broad, duplicate, equivalent and subsumed patterns; sensitivity of verdicts to molecule preparation | `smartsaudit` |
| **Screen and search** — Lipinski/Veber/Ghose/Egan/Ro3, PAINS and QED with the molecule preparation declared and recorded; substructure, similarity, fingerprints, scaffolds, R-groups | `filter` `substructure` `similar` `fingerprint` `scaffold` `rgroup` |
| **Structures** — 2D/3D generation (ETKDG + MMFF94), conformer ensembles, reaction SMILES, atom maps, SMILES augmentation | `download --gen all` `conformers` `reaction` `atommap` `augment` |
| **Databases** — one query across PubChem, ChEMBL, ChEBI and KEGG; agreement check with the *kind* of disagreement (salt form, tautomer, stereo), merged record, UniChem cross-references; name→structure concordance statistics for a whole list | `resolve` `concordance` `lookup` `batch` `download` |

## Why these features

On 9,272 ChEMBL molecules, the identical PAINS catalogue flags **77.6%, 90.5% or 88.5%** of compounds depending only on whether molecules carry implicit hydrogens, explicit hydrogens or Kekulé bonds — a choice no tool records. **16.5% of compounds change pass/fail verdict.** Of the 1,251 ChEMBL structural alerts, 358 contain a hydrogen atom and are silently dead under the default preparation, 112 are exact duplicates, and only 47 carry no flag at all.

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
