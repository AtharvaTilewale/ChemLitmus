# Compare two compound collections

**Goal:** say precisely what changed between two versions of a library — a database release, a vendor catalogue update, a generated set against its training data — in chemical terms rather than string terms.

## Why `diff` is not enough

```
-  CC(=O)Oc1ccccc1C(=O)O
+  O=C(O)c1ccccc1OC(C)=O
```

A text diff says one line was removed and one added. Chemically, nothing happened — that is aspirin written two ways. Conversely, two SMILES can differ by one character and be different drugs. Comparing compound collections needs a notion of *identity*, and the right notion depends on the question.

## Pick the level

| `--level` | Two records match when | Use when |
|---|---|---|
| `exact` | canonical isomeric SMILES identical, salts and charges included | tracking exact registrations |
| `parent` | same largest fragment after neutralisation, stereo retained | the usual default: "same active compound, ignore salt form" |
| `tautomer` | same parent, tautomer-insensitive | sources write tautomers differently |
| `nostereo` | same parent, stereo ignored | one source lacks stereo assignments |
| `skeleton` | same parent, stereo and tautomer ignored | the loosest structural match |
| `formula` | same molecular formula of the parent | isomer-level census only |

Levels nest: records that match at a stricter level also match at every looser one. See [Molecular identity](../concepts/molecular-identity.md).

## Run

```bash
chemlitmus diff chembl_34.smi chembl_35.smi --level parent --output diff.csv --json diff.json
```

```
Library Diff — chembl_34.smi → chembl_35.smi at level 'parent'
Compounds in A        2,289      3,417 valid records
Compounds in B        2,301      3,430 valid records
Added                    19      in B only
Removed                   7      in A only
Unchanged             2,241      same compound, same representation
Changed                  41      same compound at this level, different representation
Multiplicity changes      3      record count differs between A and B
Overlap (Jaccard)     0.989      shared compounds / all compounds

What changed
stereochemistry                                  25
salt, counter-ion or charge form                 12
tautomer                                          3
salt, counter-ion or charge form; stereochemistry 1
```

Four outcomes per compound:

- **added** — present in B only
- **removed** — present in A only
- **unchanged** — present in both, identical exact representation
- **changed** — present in both *at this level*, but written differently; the `change` column names the difference

"Changed" is the category a text diff cannot produce. It is where a release note gets "25 compounds had stereochemistry assigned, 12 moved to a different salt form" instead of "1,203 SMILES differ".

## The output file

`diff.csv` has one row per compound (not per record):

| Column | Meaning |
|---|---|
| `status` | `added` / `removed` / `unchanged` / `changed` |
| `level` | the identity level used |
| `key` | the identity key at that level |
| `smiles_a`, `smiles_b` | input SMILES from each file sharing this key, `|`-separated |
| `count_a`, `count_b` | how many records in each file |
| `change` | for `changed`: the reason(s) |

`unchanged` rows are omitted unless you pass `--include-unchanged`. `diff.json` carries the full `LibraryDiff` model including the summary counts.

## Interpreting multiplicity

`count_a` ≠ `count_b` for a shared compound means one file lists it more times — usually several salt forms or stereoisomers collapsing onto one key at the chosen level. The `Multiplicity changes` count summarises this. It is often the first sign that one source has started (or stopped) registering salts separately.

## Typical questions

**"Which compounds did they delete?"** — `--level parent`, read the `removed` rows.

**"Did they add stereo to existing compounds?"** — `--level nostereo`; compounds that were `removed` + `added` at `parent` become `changed: stereochemistry`.

**"Is my generated set novel?"** — `diff training.smi generated.smi --level skeleton`; `added` rows are novel even ignoring stereo and tautomer form; `unchanged`/`changed` rows are memorised, possibly with a representational twist.

**"Did the vendor silently reformat everything?"** — `--level exact` shows large churn; `--level parent` shows almost none. That gap *is* the reformatting.

## From Python

```python
from chemlitmus import diff_libraries, compute_identity

d = diff_libraries(old_smiles, new_smiles, level="parent")
print(d.n_added, d.n_removed, d.n_changed, d.jaccard)
print(d.changes_by_kind)

for e in d.entries:
    if e.status == "changed":
        print(e.smiles_a, "->", e.smiles_b, "|", e.change)

# Precompute keys once if you will diff the same collections at several levels
ka = [compute_identity(s) for s in old_smiles]
kb = [compute_identity(s) for s in new_smiles]
for lv in ("exact", "parent", "nostereo", "skeleton"):
    print(lv, diff_libraries(old_smiles, new_smiles, level=lv, keys_a=ka, keys_b=kb).n_changed)
```

Full model reference: [Identity and diff API](../reference/python-api.md#molecular-identity-and-library-comparison).
