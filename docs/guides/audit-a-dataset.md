# Audit a dataset

One command answers: what is wrong or uncertain in this dataset, how does it affect my workflow,
and what is the evidence?

The commands below run on [`examples/bioactivity.csv`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/bioactivity.csv), a 33-record table with planted problems — see [Example data](../getting-started/example-data.md) for what each record is there to show.

```bash
chemlitmus audit examples/bioactivity.csv --endpoint-column standard_value --units-column standard_units \
    --relation-column standard_relation --split-column split --date-column year -o audit_results
```

Everything runs offline. Every input row is accounted for, every original column is preserved,
and every finding carries a stable code, a severity, the affected record ids, evidence and a
suggested action.

## What it checks

| Area | Checks |
|---|---|
| Parsing | empty fields, unparseable SMILES, whitespace that RDKit would silently truncate, failed SDF records, candidate repairs |
| Standardisation | largest-fragment choice, neutralisation, optional tautomer/stereo/isotope handling — with before/after SMILES per step and whether *identity* changed, not just text |
| Structure | multiple components, more than one substantial organic component (mixtures), net charge, isotopes, unassigned or partial stereocentres, unusual elements |
| Identity | keys at all six levels, groups at the policy level, duplicate classes (`DUP_EXACT`, `DUP_PARENT`, `DUP_TAUTOMER`, `DUP_NOSTEREO`, `DUP_SKELETON`) reported separately |
| Properties | MW, logP, HBD/HBA, TPSA, rotatable bonds, heavy atoms, rings, Fsp3, with distributions |
| Alerts | configurable RDKit catalogues, with the matched atoms |
| Splits | leakage between splits when a split column is present (see [below](#leakage)) |
| Endpoints | contradictory labels or measurements when an endpoint column is present |

## Telling it about your columns

Header names are recognised automatically (`smiles`, `id`, `activity`, `split`, `units`, …). When
they are not, say so — ChemLitmus refuses to guess which of several columns holds the structures:

```bash
chemlitmus audit examples/bioactivity.csv --structure-column smiles --id-column compound_id \
    --split-column split --endpoint-column standard_value --units-column standard_units \
    --relation-column standard_relation --target-column assay_id --source-column document_id --date-column year
```

A header-less single-column file is unambiguous and needs nothing. A header-less *multi*-column
file raises a schema error unless you pass `--first-column`.

## The chemical policy

Every chemical decision is a choice; `--config` makes it explicit and hashes it into the outputs.

```bash
chemlitmus audit examples/bioactivity.csv --config conservative            # keep all components and charges, identity at 'exact', no repairs
chemlitmus audit examples/bioactivity.csv --config parent                  # largest organic fragment, neutralised, identity at 'parent' (default)
chemlitmus audit examples/bioactivity.csv --config examples/policy.json    # anything in between — this one groups at 'nostereo' and adds Brenk alerts
```

```python
from chemlitmus import ChemicalPolicy
p = ChemicalPolicy(identity_level="nostereo", preparation="explicit-h", alert_sets=["PAINS", "BRENK"],
                   fingerprint="fcfp4", similarity_threshold=0.8,
                   severity={"overrides": {"ALERT_MATCH": "ignore", "DUP_NOSTEREO": "error"}})
p.save("my_policy.json")   # p.hash identifies it exactly
```

The policy decides what counts as a duplicate, what counts as leakage, which alert sets are
applied, which fingerprint measures relatedness, whether repairs are even considered, and the
severity of every issue code.

## Outputs

| File | Contents |
|---|---|
| `records.csv` | one row per input record: every original column (prefixed `in:`), status, standardised structure, all six identity keys, structural flags, descriptors, scaffold, alerts, issue codes |
| `issues.csv` | one row per finding: code, severity, record, message, suggested action, evidence |
| `identity_groups.csv` | groups at the policy level with `differs_by` and the splits they span |
| `summary.json` | counts, denominators, gate result — the machine-readable headline |
| `audit.json` | the authoritative complete result (schema-versioned) |
| `policy.json` | the resolved policy |
| `manifest.json` | input/output checksums, tool and RDKit versions, policy hash, settings |
| `report.html` | offline report: accounting, issue cards, identity groups, descriptor histograms, leakage tables, structure depictions with alerts highlighted, searchable issue table |
| `clean.csv` + `exclusions.csv` | with `--clean`: the filtered dataset and **every** exclusion with its reason |

Nothing is dropped silently: `n_total == n_ok + n_empty + n_invalid + n_unsupported + n_error`,
and `clean.csv` + `exclusions.csv` always sum back to the input.

## Gates for CI

```bash
chemlitmus audit examples/bioactivity.csv --split-column split --fail-on-severity error --max-invalid-fraction 0.02 --no-split-overlap   # exits 3: two invalid records and a parent-level overlap
```

Exit codes: **0** pass · **1** configuration error · **3** policy violated · **4** processing
incomplete (with `--require-complete`). A failed computation is never reported as a clean pass.

```yaml
# .github/workflows/data.yml
- run: pip install chemlitmus
- run: chemlitmus audit data/train.csv --split-column split --config policy.json
        --fail-on-severity error --no-split-overlap --max-invalid-fraction 0.01
        --output-dir audit_results
- uses: actions/upload-artifact@v4
  if: always()
  with: {name: dataset-audit, path: audit_results}
```

```python
from chemlitmus import audit_dataset, GatePolicy, evaluate_gates
a = audit_dataset("data/train.csv")
gate = evaluate_gates(a, GatePolicy(fail_on_severity="error", no_split_overlap=True))
assert gate.passed, gate.violations
```

## Leakage

With a split column, the audit includes a leakage report; it also runs standalone:

```bash
chemlitmus leakage examples/bioactivity.csv --split-column split --fail-on-overlap
chemlitmus leakage train.csv test.csv valid.csv          # or one file per split; the first is the reference
```

Overlap is reported as **distinct evidence classes** — exact, parent, tautomer, nostereo,
skeleton. They nest, so **the counts are not additive**: a salt of a training compound is counted
at `parent`, `nostereo` and `skeleton`, not three times over. Beyond identity:

* **scaffold overlap** and **nearest-neighbour similarity** are *relatedness* — they say the
  evaluation is easy, not that it is invalid. Acyclic molecules are excluded from scaffold
  overlap with their count reported.
* **formula matches** are reported and explicitly **not** called leakage.
* **within-split duplicates** are counted per split.
* a **temporal** check runs when a date column is supplied, and says so when it cannot.

The report's `scope` field states what it establishes: contamination *between the collections you
supplied*. It cannot see a model's pretraining corpus, so no overlap here is not evidence that an
evaluated model never saw these compounds.

## Label conflicts

```bash
chemlitmus conflicts examples/bioactivity.csv -e standard_value --units-column standard_units \
    --relation-column standard_relation --context assay_id --tolerance 1.0
```

Records of the same compound *in the same endpoint context* are compared. Molar units are
converted with the conversion recorded (`10 uM -> 10000 nM -> p = 5.0`); the spread is measured in
log10 units (tolerance 1.0 = ten-fold). Censored values (`<`, `>`, ranges) are kept as **bounds**
and never averaged. With a units column, a measurement lacking units is flagged and excluded; with no
units column at all the endpoint is treated as dimensionless (pIC50, logP, a score) and compared on
its own scale. Nothing is resolved automatically — the suggested action is review, not "keep the most potent".

With `--source-column document_id`, each conflict group also reports **where its disagreement lives**:
`technical_spread` is the largest spread within one source (repeated entries from one experiment),
`between_source_spread` the spread of per-source medians (independent measurements that disagree).
Without source metadata the split is reported as `not determinable` rather than guessed.

## Activity cliffs and suspect labels

```bash
chemlitmus cliffs examples/bioactivity.csv -e standard_value --units-column standard_units \
    --relation-column standard_relation --context assay_id \
    -o pairs.csv --outliers outliers.csv --json cliffs.json
```

Where `conflicts` compares records of the same compound, `cliffs` compares *different* compounds
that are nearly the same — matched molecular pairs (one site changed, including H → substituent)
and fingerprint neighbours at Tanimoto ≥ 0.9 — and reports every pair whose labels are **proven**
to differ by at least the threshold (default 1.0 log unit = ten-fold). Censored values are bounds
here too: a pair that cannot be decided is `undetermined`, never a cliff, and the cliff fraction is
quoted over decided pairs only.

Compounds whose every neighbour disagrees with them while the neighbours agree with each other are
listed as **label outliers** (`LABEL_OUTLIER`, warning). On 6,000 ChEMBL IC50 records for
acetylcholinesterase this surfaces, among others, a compound reported at pIC50 11.2 whose three
one-atom neighbours sit at 7.8–8.7 — the signature of a unit error. The command does not decide:
a genuine cliff produces the same pattern, and the suggested action is to check the source record.

The `transformations` table groups matched pairs by the fragment exchange (`C[*:1]>>[H][*:1]`,
…) with pair counts by verdict and the mean signed difference, which separates systematic SAR from
isolated disagreements. See [Activity cliffs and label outliers](../concepts/activity-cliffs.md).

## Then: splits and generated molecules

```bash
chemlitmus split examples/bioactivity.csv -s scaffold -f train=0.8,test=0.2 --seed 0 -o split.csv
chemlitmus generated examples/generated.smi -r examples/train.smi --constraints 'mw=0:500,logp=-1:5'
```

See [Build splits and evaluate generated molecules](splits-and-generation.md).
