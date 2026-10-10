# Example data

All example files live in the repository's [`examples/`](https://github.com/AtharvaTilewale/ChemLitmus/tree/main/examples) directory. Every command in these docs that takes an input file is written against one of them, so each can be run as shown.

Small, hand-built input files that every command in the documentation runs on. Each file is
deliberately imperfect — the defects are the point — and every row that is "wrong" is wrong on
purpose, so the documented output can be reproduced exactly.

```bash
git clone https://github.com/AtharvaTilewale/ChemLitmus.git
cd ChemLitmus
chemlitmus audit examples/bioactivity.csv --endpoint-column standard_value --units-column standard_units
```

Individual files can also be fetched directly, e.g.
`curl -O https://raw.githubusercontent.com/AtharvaTilewale/ChemLitmus/main/examples/bioactivity.csv`.

| File | Rows | What it is | Used by |
|---|---|---|---|
| [`bioactivity.csv`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/bioactivity.csv) | 33 | A bioactivity table in ChEMBL-style columns (`compound_id, smiles, standard_relation, standard_value, standard_units, pchembl_value, assay_id, document_id, year, split, activity_class`). Two SAR series, salt / stereo / tautomer forms of the same compound spread across `train` and `test`, one compound measured 1,600-fold apart in two documents, a censored `>` value, a value with no units, a record with no split label, a mixture, isotopes, partial stereo, a permanent charge, a platinum complex, a PAINS quinone, and three records that do not parse (unbalanced ring, internal whitespace, empty). | `audit`, `leakage`, `conflicts`, `cliffs`, `split` |
| [`library.csv`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/library.csv) | 24 | A vendor library (`id, smiles, name, vendor`): drugs written as salts, Kekulé forms, enantiomers, tautomer pairs and exact duplicates, plus a PAINS catechol and quinone, a very lipophilic acid, a permanently charged compound, and three broken records. | `validate`, `diagnose`, `standardize`, `identity`, `filter`, `similar`, `substructure`, `fingerprint`, `scaffold`, `tautomers` |
| [`library.smi`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/library.smi) | 23 | The same library as SMILES-plus-id lines. | `substructure --prep explicit-h`, `download --file`, `validate --file` |
| [`library_v2.csv`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/library_v2.csv) | 21 | A later release of the same library: two compounds removed, two added, ibuprofen switched to its sodium salt, the broken rows cleaned up, aspirin now listed by a second vendor. | `diff` |
| [`demo.csv`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/demo.csv), [`demo_v2.csv`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/demo_v2.csv) | 9, 5 | The nine-record file used in the [Quickstart](quickstart.md) and its revised version. | Quickstart |
| [`alerts.csv`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/alerts.csv) | 18 | A structural-alert catalogue (`name, smarts, rule_set`) with real alerts and planted defects: an unparseable pattern, an exact duplicate, an aromatic/Kekulé equivalent pair, a pattern that can never match (`[#118]`), one that needs explicit hydrogens (`c[H]`), an over-broad `[#6]`, a recursive SMARTS, and several patterns provably contained in others. | `smartsaudit`, `smartsproof`, `smartsdiff` |
| [`alerts_v2.csv`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/alerts_v2.csv) | 13 | A revised catalogue: six patterns removed (the over-broad `[#6]` and `c1ccccc1`, the duplicate, the Kekulé twin, the broken pattern and `[#118]`), two broadened (`aldehyde`, `michael_acceptor`), two rewritten without changing their hits on the reference panel (`phenol`, `thiol`), one added (`epoxide`). | `smartsdiff` |
| [`alerts.smarts`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/alerts.smarts) | 5 | The five-line catalogue from the Quickstart, one SMARTS and name per line. | `smartsaudit` |
| [`train.smi`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/train.smi) | 15 | A "training set" of known drugs. | `generated --reference`, `leakage` |
| [`generated.smi`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/generated.smi) | 20 | Raw output of an imaginary generator: two unparseable strings, an exact duplicate, two tautomers of one compound, two copies of a training compound (one as a salt), a PAINS quinone, a 46-carbon alkane that fails the MW constraint, and otherwise novel molecules. | `generated` |
| [`analogues.smi`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/analogues.smi) | 6 | Mono-substituted benzenes for R-group decomposition against `c1ccccc1[*:1]`. | `rgroup` |
| [`names.txt`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/names.txt) | 8 | Drug names for the network-backed lookups. | `resolve --file`, `concordance --file` |
| [`smiles.txt`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/smiles.txt) | 5 | SMILES for batch metadata retrieval. | `batch` |
| [`policy.json`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/policy.json) | — | A chemical policy that groups at `nostereo`, applies PAINS + Brenk alerts, lowers the similarity threshold to 0.8, promotes `DUP_NOSTEREO` to a warning and ignores `ISOTOPE`. | `--config examples/policy.json` on any dataset command |

## What the dataset commands find in `bioactivity.csv`

| Command | Finding | Records |
|---|---|---|
| `audit` | 33 records → 30 parsed, 1 empty, 2 invalid; 28 distinct compounds at `parent` | — |
| `audit` / `leakage` | aspirin free acid (train) and sodium salt (test) share a parent; caffeine appears in both splits exactly; racemic vs (S)-ibuprofen share a `nostereo` key; 2-hydroxypyridine / 2-pyridone share a `tautomer` key | CPD-015/016, CPD-017/018, CPD-019/020, CPD-021/022 |
| `conflicts` | caffeine measured at 80 µM in DOC-3 and 50 nM in DOC-4 — a 3.2 log-unit spread between sources | CPD-017, CPD-018 |
| `cliffs` | the nitro benzanilide at 4 nM against nine one-site neighbours at 2–15 µM (and one `> 100 µM`): a **label outlier** with 8 of 9 neighbours agreeing; the amino- and methylamino-sulfonamides at 0.3–0.5 µM against their H and Cl analogues at 40–50 µM: a **genuine cliff** supported by two compounds | CPD-007; CPD-013/014 vs CPD-011/012 |
| `audit` | missing units, censored relation, mixture, isotopes, partial stereo, permanent charge, unusual element, PAINS match, no split label | CPD-030, CPD-008, CPD-024, CPD-025, CPD-023, CPD-027, CPD-028, CPD-026, CPD-029 |

The example files are fictional. Activity values were chosen to illustrate behaviour and are not
measurements of the named compounds.
