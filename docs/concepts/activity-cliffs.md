# Activity cliffs and label outliers

`conflicts` compares records of the **same** compound. `cliffs` compares **different** compounds
that are nearly the same, and asks whether their labels agree. The two questions are complementary:
a dataset can be free of duplicate conflicts and still contain values that no neighbouring
measurement supports.

## Why this matters

A large label difference across a small structural change is an *activity cliff*. Cliffs are real
chemistry — a methyl that fills a pocket, a halogen that breaks a hydrogen bond — and the pairs that
exhibit them carry more structure–activity information than any other pairs in the set. They are
also exactly what a wrong value looks like: a unit mistake, a transcription error or a mis-assigned
structure turns a compound into an apparent cliff against every neighbour it has.

For machine learning the consequences run both ways. Models fit on cliff-rich series look worse
than they are if the cliffs are real (the test asks a genuinely hard question) and better than
they are if the cliffs are errors that the model learned to reproduce. Benchmarks such as
MoleculeACE exist because neither the error rate nor the difficulty of a set can be read from
summary statistics; both require looking at pairs.

`cliffs` reports every pair with its evidence and never decides which explanation holds. It
additionally flags *label outliers* — the pattern a wrong value produces most often — as review
candidates, not as errors.

## How compounds are paired

Two methods, each reported by name on every pair (`relationship`: `mmp`, `similarity` or
`mmp+similarity`):

**Matched molecular pairs.** Each compound is cut once at every acyclic single bond (RDKit's
`rdMMPA` fragmentation) and, additionally, at every hydrogen, so that benzene / toluene and
N–H / N–methyl are pairs. The two halves are a *core* (with the attachment point written `[*:1]`)
and a *varied fragment*. Following Hussain & Rea (2010), the varied fragment may be at most
`--max-r-atoms` heavy atoms (default 13) and never larger than the core. Two compounds that share a
core with different varied fragments are a pair; when several cores are shared, the largest
(smallest change) is reported, as `transformation` `R_a>>R_b` with `n_changed_atoms`.

**Similarity.** Tanimoto on the policy fingerprint at or above `--similarity-threshold` (default
0.9). This catches changes at two sites, ring-size changes and other single edits that the one-cut
rule misses; it is also reported on MMP pairs for information.

Pairs are formed only within the same endpoint context (`--context target,assay`), and only between
compounds that are **different at the policy identity level** — records of the same compound are
merged into one value first.

## How labels are compared

Classification labels agree or they do not.

Quantitative values are compared on one scale: log10 molar (pX) when every usable value has a molar
unit, otherwise the column's own units — or, when there is no units column at all, the raw numbers
as a dimensionless quantity (pIC50, logP, a score). The threshold (default 1.0) is a difference on
that scale: ten-fold for molar data.

Censored values are **bounds**, exactly as in `conflicts`. `> 10 µM` is the interval pX < 5;
`< 1 nM` is pX > 9; a range is the interval between its ends. Each pair then has a *proven* minimum
difference (the gap between the two intervals, 0 when they overlap) and a maximum compatible
difference (open when either bound is open). The verdict follows from those two numbers and nothing
else:

| Verdict | Condition |
|---|---|
| `cliff` | the minimum proven difference is at or above the threshold |
| `consistent` | the maximum compatible difference is below the threshold |
| `undetermined` | neither can be shown — typically one or both values censored |

Two compounds both reported as `> 50 µM` are *undetermined*, not consistent: each could lie anywhere
below pX 4.3. The `cliff_fraction` is computed over decided pairs only and its denominator is stated.

A compound whose own records disagree beyond the threshold is excluded (`n_internal_conflict_excluded`)
and left to `conflicts`; one with several agreeing records contributes their median, labelled as such.

## Label outliers

A compound is a label outlier when

1. every near neighbour it has is a proven cliff against it, and
2. at least `--min-neighbours` of those neighbours (default 2) are proven consistent with one
   another.

Condition 2 is a clique over exact or tightly bounded values. Neighbours that disagree among
themselves — a genuine cliff inside the series, or a censored value that cannot be decided — do
not block the call, but they do not vouch for anything either: the report gives `n_neighbours`
and `n_agreeing` so the strength of the call is visible (9 neighbours of which 8 agree is a
different statement from 2 of 2).

The call is an issue of severity `warning` (`LABEL_OUTLIER`) because the pattern is the one a
wrong value produces most often — but a genuine cliff at the edge of a series produces it too. The
suggested action is to check the source record, not to change it.

## Transformations

Matched pairs are also grouped by the fragment exchange that relates them, written in a canonical
orientation (`R1>>R2` with `R1` sorting first), with pair counts by verdict and the mean *signed*
difference over pairs with exact values. A transformation whose pairs are mostly cliffs with a
consistent sign is SAR; a transformation that produces a single cliff among many consistent pairs
points back at that one pair.

## What this is not

* **Not a model.** No activity is predicted; a neighbour's value is evidence, not a prediction.
* **Not a verdict on correctness.** Both members of a cliff may be right.
* **Not a cleaning step.** Nothing is removed or altered; the outputs are a pair table, an outlier
  table and issue codes (`ACTIVITY_CLIFF` info, `LABEL_OUTLIER` warning) whose severities the policy
  may override.
* **Not exhaustive.** Single-cut pairs and a fingerprint threshold define *near*; compounds that
  differ in more than that are outside the analysis, and a set with no pairs has not been shown to be
  cliff-free.

## References

* Hussain J, Rea C. *J Chem Inf Model* 2010, 50, 339–348 — matched molecular pair conventions.
* Stumpfe D, Bajorath J. *J Med Chem* 2012, 55, 2932–2942 — activity cliffs.
* van Tilborg D, Alenicheva A, Grisoni F. *J Chem Inf Model* 2022, 62, 5938–5951 — MoleculeACE.
