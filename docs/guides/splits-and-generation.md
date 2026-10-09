# Build splits and evaluate generated molecules

## Group-aware splits

```bash
chemlitmus split dataset.csv --strategy identity --fractions train=0.8,test=0.2 --seed 0 -o split.csv
```

| Strategy | Group = |
|---|---|
| `random` | the record itself |
| `identity` | identity key at the policy level (a salt and its free acid stay together at `parent`) |
| `scaffold` | Murcko scaffold; **each acyclic molecule is its own group**, never pooled into one giant group |
| `temporal` | the date value (needs `--date-column`); earlier groups fill the earlier splits |
| `source` | the source/document value (needs `--source-column`) |

**Groups are never divided.** That makes achieved fractions approximate, so the report states the
achieved fractions, the group-size distribution (including the largest group as a share of the
data), the endpoint balance per split, and every constraint that could not be satisfied:

```
  constraint: largest group holds 8 records (88.9%), more than the smallest requested fraction (50.0%);
              that split is necessarily larger than requested
  constraint: split 'test': achieved 11.1% vs requested 50.0% (indivisible groups)
```

The split is then **verified with the same leakage implementation** the audit uses, so the
report shows the overlap it achieved rather than the overlap it intended. Seeds reproduce
`random`, `identity` and `scaffold` assignments exactly.

A scaffold split is an evaluation design, not a guarantee of deployment generalisation.

## Evaluating generated molecules

```bash
chemlitmus generated samples.smi --reference train.smi --constraints 'mw=0:500,logp=-1:5' --repair
```

Pass the raw outputs — do not pre-filter. Every metric states its denominator:

* **validity** = parseable / **all generated attempts** (empty strings included, and counted);
* **uniqueness** = distinct keys / **valid outputs**, reported at each identity level (they nest,
  so they are not independent);
* **novelty** = outputs absent from the reference sets / **unique valid outputs**, and claimed
  only against the collections actually supplied — never as novelty in general.

Also reported: scaffold diversity, nearest-reference Tanimoto (fingerprint named), descriptor
summaries, alert prevalence with its denominator, and the share of valid outputs satisfying each
explicit constraint.

Novelty depends on the identity policy, which is the point of making it explicit:

```python
from chemlitmus import ChemicalPolicy, evaluate_generated
gen = ["CC(=O)Oc1ccccc1C(=O)[O-].[Na+]", "C[C@@H](N)C(=O)O"]
ref = {"train": ["CC(=O)Oc1ccccc1C(=O)O", "C[C@H](N)C(=O)O"]}
evaluate_generated(gen, ref, ChemicalPolicy(identity_level="exact")).novelty      # 1.0  — both strings differ
evaluate_generated(gen, ref, ChemicalPolicy(identity_level="parent")).novelty     # 0.5  — the salt is not a new compound
evaluate_generated(gen, ref, ChemicalPolicy(identity_level="nostereo")).novelty   # 0.0  — nor is the enantiomer
```

**Repaired candidates are a separate population.** With `--repair`, mechanical repairs of invalid
outputs are evaluated in their own report; they are never folded into the model's validity.

When nothing valid was generated, metrics that have no denominator are listed under `undefined`
rather than reported as zero:

```python
evaluate_generated(["bad((", ""]).undefined
# ['uniqueness', 'novelty', 'scaffold_diversity', 'nearest_neighbour_similarity']
```

Descriptor compliance and alert counts describe structure. They are not evidence of potency,
safety, ADMET behaviour or synthesisability.
