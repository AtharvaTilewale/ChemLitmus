# Static SMARTS containment proofs

`smartsaudit` finds patterns that are *empirically* redundant: on the reference library, one
pattern's hit set contains another's. That is evidence, not proof — a different library could
separate them. `smartsproof` proves containment **for every molecule**, from the patterns alone.

## The claim

**P ⊆ Q** ("Q contains P", "P is redundant given Q") means: every molecule that matches P also
matches Q. If both are in one catalogue, P adds no compound to the flagged set.

## How a proof is built

1. **Atom expressions** are read from RDKit's query tree (the output of `DescribeQuery()`), so
   the prover reasons about exactly what RDKit will match — not about SMARTS text.
2. Each expression is evaluated over a **finite universe of atom states**: element ×
   aromaticity × formal charge × total H count × explicit degree × implicit H count × ring count
   × ring-bond count. The universe is built per run from the values the patterns actually
   mention, plus one unreferenced "other" value in every dimension, and pruned by the physical
   constraints that always hold (ring bonds ≤ degree; in a ring ⇔ ≥ 2 ring bonds; aromatic ⇒ in
   a ring; total H ≤ implicit H + degree). Over that universe, "A implies B" is decided exactly:
   `bits(A) & ~bits(B) == 0`. Bond expressions are handled the same way over bond order × ring
   membership.
3. **P ⊆ Q is proven by a witness**: an injective map from the atoms of Q into the atoms of P
   such that every P atom implies the Q atom mapped onto it, and every Q bond lands on a P bond
   whose expression implies it. This is the Chandra–Merlin homomorphism criterion for query
   containment, specialised to RDKit's injective, non-induced substructure semantics. The witness
   is reported so the claim can be checked by hand.

```bash
chemlitmus smartsproof --subsumes 'c1ccccc1[OH]' 'c[OH]'
# A ⊆ B: proven   Witness: B[0]→A[5], B[1]→A[6]
chemlitmus smartsproof --subsumes 'CC' 'CCO'
# A ⊆ B: no witness  — not a refutation
```

## What the answers mean

| Status | Meaning |
|---|---|
| `proven` | A witness exists. Containment holds for every molecule. |
| `no witness` | **Undecided.** The method is sound but not complete; a containment may hold without a witness of this form. Never read this as "not contained". |
| `budget exceeded` | The search stopped early (container larger than `--max-container-atoms`, or more than `--step-budget` backtracking steps). Undecided. |
| `not analysable` | The pattern uses a primitive outside the universe (below). |
| `unsatisfiable` | An atom expression matches no atom state — `[C;N]`, `[c;!a]`, `[R0;x2]`. The pattern can never fire. |

Every count of proven redundancy is therefore a **lower bound**.

## Limits

* **Outside the universe**: recursive SMARTS `$(...)`, isotopes (`[2H]`, `[13C]`), total
  valence `v<n>`, ring size `r<n>`, chirality, and disconnected / component-grouped patterns.
  On the ChEMBL alert collection this excludes 216 of 1,251 patterns (182 of them recursive).
* **Default match semantics.** RDKit ignores chirality and bond direction unless asked; so does
  the prover. Containment under `useChirality=True` is not what is proven.
* **Satisfiability is local.** `[N+]#[C-]` passes: each atom and bond is individually
  realisable. Global impossibilities (valence, ring geometry) are not detected.
* **Redundancy is not always a defect.** A pattern duplicated across *published* sets carries
  provenance; the tool reports and lets you decide.
* **Soundness depends on the universe covering every realisable state.** Elements, charges,
  degrees and H counts that no pattern mentions fall into the "other" value; a state that a
  pattern does mention is enumerated. Hypervalent or exotic states are represented by the
  degree and charge values present, which is why `C(=O)[OH] ⊆ [CX3](=O)[OX2H1]` is *not*
  proven: nothing in the left pattern bounds the carbon's degree.

## Results on the ChEMBL structural-alert collection

Running `chemlitmus smartsproof` on the 1,251 alerts (eight published sets) takes about ten
seconds: 557,460 atom states, 552,437 candidate pairs.

| Outcome | Patterns |
|---|---|
| Proven redundant (another alert provably contains it) | **395 (31.6%)** |
| of which proven equivalent | 187 |
| Undecided (no witness) | 640 |
| Not analysable | 216 (182 recursive, 26 disconnected, 5 valence, 5 isotope, 2 ring size) |
| Unsatisfiable | 0 |

Of the 395, a proven container lies in a *different* published set for 386 and in the same set
for 49 (some have both): this is the same chemistry re-encoded by several groups, far more than
duplication within one set.
The `smartsaudit --checks all,proof` run places the proofs next to the empirical verdicts:
`proven_subsumed_by` versus `subsumed_by`, so you can see which empirical claims are proven,
which are undecided, and which empirical redundancies the prover cannot speak to.
