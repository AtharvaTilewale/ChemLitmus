# Molecular identity

"Are these two records the same compound?" has no single answer. A sodium salt and its free acid are the same active ingredient and different substances. Two tautomers are one compound in solution and two in a crystal. A racemate and one enantiomer share a connection table and differ in biology. The right answer depends on what you are doing — so ChemLitmus does not pick one. It computes identity at six nested levels and lets you choose.

## The six levels

| Level | What is kept | What is ignored | Computed as |
|---|---|---|---|
| `exact` | everything in the input | nothing | canonical isomeric SMILES |
| `parent` | largest organic fragment, neutralised, with stereo | counter-ions, charge state | canonical SMILES of the parent |
| `tautomer` | parent, with stereo | + tautomeric form | RDKit `RegistrationHash` `TAUTOMER_HASH` of the parent |
| `nostereo` | parent, with tautomer | + stereochemistry | `NO_STEREO_SMILES` of the parent |
| `skeleton` | parent | + stereo + tautomer | `NO_STEREO_TAUTOMER_HASH` of the parent |
| `formula` | element counts of the parent | everything else | `FORMULA` of the parent |

The levels nest: two records identical at a stricter level are identical at every looser one. `tautomer` and `nostereo` are siblings — each coarsens `parent` along one axis, and `skeleton` coarsens along both.

```
exact ⊂ parent ⊂ { tautomer , nostereo } ⊂ skeleton ⊂ formula
```

The parent is produced by RDKit's `LargestFragmentChooser` followed by `Uncharger` — the same two steps `standardize` uses for `fragment` and `neutralize`. The hashed levels use RDKit's `RegistrationHash` module, which is the machinery compound-registration systems use to decide whether a new submission is already in the database. ChemLitmus exposes that decision at every resolution instead of fixing one.

## Choosing a level

| Question | Level |
|---|---|
| Is this exact registration already present? | `exact` |
| Same active compound, ignoring how it was salted? (the usual default) | `parent` |
| Did two sources draw the same compound as different tautomers? | `tautomer` |
| One source has stereo assigned and the other does not — same compound? | `nostereo` |
| Loosest structural match: same 2D skeleton regardless of stereo or tautomer | `skeleton` |
| Isomer census only | `formula` |

For deduplicating a screening library, `parent` is almost always right. For deciding whether a generated molecule is novel against a training set, `skeleton` is the stricter test of novelty (it refuses to call a stereoisomer or tautomer of a training compound "new"). For reconciling databases that disagree about tautomer drawing conventions, `tautomer`.

## Reading the identity report

```
Level      Distinct compounds  Collapsed records
exact                   12,480                  0
parent                  11,902                578
tautomer                11,715                765
nostereo                11,340              1,140
skeleton                11,201              1,279
```

Each step down the table is a question answered: 578 records are salt or charge variants of another record; a further 187 are tautomer variants; stereo variation accounts for 562; and so on. The groups table then shows each cluster with a `differs_by` list, so you can see *which* axis each duplicate lies on and act accordingly — dedupe salts but keep enantiomers, for instance.

## Describing a pairwise difference

`describe_difference(a, b)` reports the strictest level at which two molecules agree, translated into words:

| Shared level | Description |
|---|---|
| `exact` | identical |
| `parent` | salt, counter-ion or charge form |
| `tautomer` | tautomer |
| `nostereo` | stereochemistry |
| `skeleton` | stereochemistry and tautomer |
| `formula` | constitution (same formula only) |
| none | different compounds |

This is what `diff` writes into its `change` column.

## Limits

- **Parent selection is a heuristic.** "Largest organic fragment" is RDKit's rule; for a 1:1 co-crystal of two drug-sized molecules it will keep one and discard the other. Check `n_fragments` on records where this matters.
- **Neutralisation is chemically reasonable, not exhaustive.** Quaternary ammoniums stay charged; zwitterions are handled by RDKit's `Uncharger` rules, which do not cover every case.
- **Tautomer hashing follows RDKit's tautomer model.** It is designed to be stable, not to enumerate every chemically accessible tautomer. Two forms RDKit does not consider tautomers will have different `tautomer` keys even if a chemist would call them the same.
- **Stereo removal is total at `nostereo`/`skeleton`.** Partial stereo specification (one centre assigned, one not) collapses onto the fully unspecified form.
- **Validity is a prerequisite.** Records that fail `validate` — including those with internal whitespace — are excluded from grouping and counted as invalid in the report.
