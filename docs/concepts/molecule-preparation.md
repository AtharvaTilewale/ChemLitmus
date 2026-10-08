# Molecule preparation

Substructure matching compares a query pattern to a molecular graph. What that graph *looks like* — whether hydrogens are atoms or annotations, whether a benzene ring has alternating double bonds or six aromatic bonds — is a choice the toolkit makes before matching begins. Different choices give different matches. The choice is almost never recorded.

## The three preparations

| `--prep` | Hydrogens | Aromatic rings | Produced by |
|---|---|---|---|
| `implicit-h` | implicit (atom property, not graph nodes) | aromatic bonds, lower-case atoms | RDKit `MolFromSmiles` default |
| `explicit-h` | explicit graph nodes | aromatic | `Chem.AddHs(mol)` |
| `kekule` | implicit | alternating single/double, aromatic flags cleared | `Chem.Kekulize(mol, clearAromaticFlags=True)` |

## Why it changes the answer

**Hydrogen primitives.** A SMARTS atom written `[H]` or `[#1]` is a hydrogen *atom*. Under `implicit-h` there are no hydrogen atoms in the graph, so the pattern can never match. The PAINS set as distributed by ChEMBL has 358 such patterns — 29% of the catalogue — and every one is silently dead unless the molecules carry explicit hydrogens. The correct implicit-H way to say "carbon with two hydrogens" is `[CH2]`, an H-count query on the carbon.

**Kekulé versus aromatic.** The pattern `C1=CC=CC=C1` describes a six-membered ring of aliphatic carbons with alternating double bonds. An aromatic-perceived benzene has six *aromatic* carbons with aromatic bonds. They do not match. `c1ccccc1` matches aromatic benzene; `C1=CC=CC=C1` matches kekulised benzene; neither matches the other. A catalogue written against one representation fails against the other.

**Bond-order queries inside rings.** Patterns such as `[#6]=[#8]` inside an aromatic ring behave differently depending on whether the ring was kekulised; the exocyclic carbonyl of a pyridone is aromatic-adjacent under one preparation and a plain double bond under another.

## How much it matters

ChemLitmus measured this on the 1,251 ChEMBL structural alerts against its 9,272-molecule reference set:

| Preparation | Compounds flagged | Alerts that fire | Compounds whose verdict differs from `implicit-h` |
|---|---|---|---|
| `implicit-h` | 7,195 (77.6%) | 523 | — |
| `explicit-h` | 8,392 (90.5%) | 519 | 1,527 (16.5%) |
| `kekule` | 8,207 (88.5%) | 422 | 1,538 (16.6%) |

One compound in six changes pass/fail depending on a choice that is never written down. Two groups applying "the PAINS filter" to the same library will disagree on roughly a sixth of it, and neither will know.

## What ChemLitmus does about it

Three things:

1. **`filter` and `substructure` take `--prep`** and refuse unknown values. The default is `implicit-h`, matching RDKit, so existing behaviour is unchanged — but now it is explicit.
2. **Every output row records the preparation** (`preparation` column in CSV, `.preparation` on `FilterResult` and `SubstructureHit`). A screen whose preparation is in the file can be reproduced; one whose preparation is in someone's head cannot.
3. **`smartsaudit` measures sensitivity** for any catalogue — per pattern (does its hit count change?) and per catalogue (how many compounds flip verdict?) — and flags hydrogen-bearing patterns as `needs-explicit-h` so they are never silently dead again.

## Choosing a preparation

- If the catalogue's documentation says — use that.
- If it does not, run `smartsaudit`. If many patterns are `needs-explicit-h`, the catalogue was written for `explicit-h`. If many are `fires-only-with-kekule`, it was written against kekulised structures.
- If nothing points either way, use `implicit-h`, because it is what everyone else is unknowingly using, and **record it**.

The point is not that one preparation is right. The point is that an unrecorded preparation makes the result unreproducible, and a recorded one does not.
