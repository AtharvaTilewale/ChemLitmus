# Quickstart

Ten minutes, one small file, every major capability. Make a working directory and save this as `demo.csv`:

```csv
smiles,name
CC(=O)Oc1ccccc1C(=O)O,aspirin
CC(=O)Oc1ccccc1C(=O)[O-].[Na+],aspirin sodium
C[C@H](N)C(=O)O,L-alanine
C[C@@H](N)C(=O)O,D-alanine
Oc1ccccn1,2-hydroxypyridine
O=c1cccc[nH]1,2-pyridone
C1CC(C,broken
c1cncc1,pyrrole-typo
CC O,spaced
```

Nine records. Three of them are broken in different ways, and the other six are really three compounds written two ways each.

## 1. Validate

```bash
chemlitmus validate --file demo.csv
```

Six valid, three invalid. `validate` is the fast gate — exit code 0/2 in single mode makes it a one-line guard in a shell script:

```bash
chemlitmus validate "CC O" --quiet || echo "rejected"
```

## 2. Diagnose the failures

`validate` says *that* a SMILES failed. `diagnose` says *where* and *why*, and tries to fix it:

```bash
chemlitmus diagnose --file demo.csv
```

```
  C1CC(C
   ^  ^
  [parentheses] '(' opens a branch that is never closed  (position 4)
  [rings] Ring closure 1 is opened but never closed  (position 1)
  Repaired: CCCC

  c1cncc1
  ^
  [aromaticity] Can't kekulize mol.  Unkekulized atoms: 0 1 2 3 4
      → A pyrrole-type nitrogen must carry its hydrogen: write it as [nH]
  Repaired: c1cc[nH]c1

  CC O
    ^
  [characters] Whitespace splits the record: RDKit silently parses only 'CC' …
  Repaired: CCO
```

The third case matters more than it looks: `CC O` is **not** rejected by RDKit — it parses as ethane and quietly discards the rest. `diagnose` catches this; a plain validity check does not. Repairs are mechanical, not chemical; treat them as candidates to review.

## 3. How many compounds is this, really?

```bash
chemlitmus identity --file demo.csv --level parent
```

```
Level      Distinct compounds  Collapsed records
exact                       6                  0
parent <                    5                  1
tautomer                    4                  2
nostereo                    4                  2
skeleton                    3                  3
```

Six valid records; **three compounds**. At `parent` level the sodium salt collapses onto aspirin. At `nostereo` the two alanines merge. At `tautomer` the two pyridines merge. At `skeleton` all three distinctions are ignored. The groups table tells you *what* varies within each group — salt form, stereochemistry, tautomer — so you can decide which level your task needs. See [Molecular identity](../concepts/molecular-identity.md).

## 4. Standardise

```bash
chemlitmus standardize "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]" --show-diff
```

Strips the counter-ion, neutralises the carboxylate, canonicalises the tautomer, and shows each step. Batch mode writes a CSV with one row per input.

## 5. Compare two versions of a library

Save this as `demo_v2.csv` — D-alanine and 2-hydroxypyridine are gone, the broken rows are gone, ethylamine is new:

```csv
smiles,name
CC(=O)Oc1ccccc1C(=O)O,aspirin
CC(=O)Oc1ccccc1C(=O)[O-].[Na+],aspirin sodium
C[C@H](N)C(=O)O,L-alanine
O=c1cccc[nH]1,2-pyridone
CCN,ethylamine
```

```bash
chemlitmus diff demo.csv demo_v2.csv --level parent
```

At `parent` level: **1 added**, **3 removed**, **3 unchanged**. The three removals are D-alanine, 2-hydroxypyridine, and — because the salt is already folded in — nothing else surprising. Now change the level:

```bash
chemlitmus diff demo.csv demo_v2.csv --level nostereo    # removed 2, changed 1
chemlitmus diff demo.csv demo_v2.csv --level skeleton    # removed 1, changed 2
```

At `nostereo`, D-alanine is no longer "removed" — it is the same compound as L-alanine at that resolution, so the alanine entry becomes **changed: stereochemistry**. At `skeleton`, the dropped hydroxypyridine likewise becomes **changed: tautomer** of the pyridone that survived. Same two files, three honest answers, depending on what you mean by "the same compound". A text `diff` cannot make that distinction at all.

## 6. Screen — and say how

```bash
chemlitmus filter --file demo.csv --rules lipinski,pains --prep explicit-h --output screened.csv
```

The `--prep` flag is new and worth understanding. Structural-alert matching gives *different answers* depending on whether molecules carry explicit hydrogens or Kekulé bonds — on a 9,272-molecule reference set the identical PAINS catalogue flags 78% of compounds one way and 91% another. `filter` and `substructure` therefore require the preparation to be declared and write it into every output row, so the screen can be reproduced. See [Molecule preparation](../concepts/molecule-preparation.md).

## 7. Audit the filter itself

Nobody checks structural-alert catalogues. ChemLitmus does. Save a few patterns as `alerts.smarts`:

```
c1ccccc1   benzene
[#6]       any-carbon
c1ccccc1   benzene-again
[#118]     oganesson
c[H]       aromatic-CH
```

```bash
chemlitmus smartsaudit alerts.smarts
```

The audit reports that `benzene-again` is an exact duplicate, that `benzene` is subsumed by `any-carbon` (every benzene contains a carbon, so the stricter rule can never change a rejection decision), that `oganesson` contains an atom no real molecule has, and that `aromatic-CH` is dead under the default preparation but fires with explicit hydrogens. Point it at a real catalogue — the ChEMBL structural alerts, your in-house filters — and the picture is the same, at scale. See [Auditing alert sets](../guides/audit-alert-sets.md).

## 8. From Python

Everything above is a function returning a typed model:

```python
from chemlitmus import diagnose_smiles, group_by_identity, diff_libraries, audit_smarts

d = diagnose_smiles("CC O")
print(d.is_valid, d.problems[0].message, d.repaired_smiles)

rep = group_by_identity(["C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"], level="nostereo")
print(rep.n_groups, rep.groups[0].differs_by)      # 1  ['stereochemistry']
```

## Next steps

- [Guides](../guides/clean-a-library.md) — end-to-end workflows
- [CLI reference](../reference/cli.md) — every option
- [Configuration](configuration.md) — cache location, PubChem rate limits, logging
