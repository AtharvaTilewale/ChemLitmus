# Audit a structural-alert set

**Goal:** before using a SMARTS catalogue — PAINS, Brenk, Glaxo, an in-house filter — to reject compounds, find out which of its rules are broken, dead, redundant or unstable, and under which molecule preparation it was meant to be run.

## Why this is necessary

Structural-alert catalogues are redistributed as plain text and applied with whatever molecule preparation a toolkit defaults to. Nobody validates them. When ChemLitmus audited the 1,251 ChEMBL structural alerts (PAINS, BMS, SureChEMBL, MLSMR, Dundee, Inpharmatica, LINT, Glaxo) against its 9,272-molecule reference set:

- Taken together, the eight sets flagged **77.6%, 90.5% or 88.5%** of compounds depending only on whether molecules had implicit hydrogens, explicit hydrogens or Kekulé bonds; **16.5% of compounds changed pass/fail verdict** between implicit and explicit H. *Per set the spread differs enormously* — SureChEMBL flips 42.5% of compounds, PAINS 1.3% — so never quote a union figure for one named set (see [Molecule preparation](../concepts/molecule-preparation.md#how-much-it-matters)).
- **358 alerts (29%) contain a hydrogen atom** and therefore match nothing under RDKit's default preparation. They are silently dead unless you add hydrogens.
- **112 are byte-identical duplicates** of another alert; **354 are empirically subsumed** by another (every molecule they flag is already flagged by a broader rule).
- Only **47 of 1,251 (3.8%)** carried no flag at all.

All of these are measurements on *this* catalogue file against *this* reference panel; the recipe is [`benchmarks/alert_preparation_sensitivity.py`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/benchmarks/alert_preparation_sensitivity.py), which prints the catalogue checksum alongside every figure.

None of this is visible from the SMARTS text. It is only visible by measurement.

## Input format

A CSV/TSV/XLSX with a column named `smarts` (or `pattern`, `query`), optionally `description`/`name` and `rule_set_name`/`rule_set`:

```csv
rule_set_name,description,smarts
PAINS,azo_A(324),c:1:c:c(:c:c:c:1-[#7]=[#7]-c:2:c:c:c:c:c:2)-[#8]-[#1]
Glaxo,R1 Reactive alkyl halides,"[Br,Cl,I][CX4;CH,CH2]"
```

Or a plain text file, one SMARTS per line with an optional name after whitespace; `#` starts a comment.

The ChEMBL alert table as redistributed by `rd_filters` works directly: `alert_collection.csv` has exactly these columns.

## Run the audit

```bash
chemlitmus smartsaudit examples/alerts.csv --output audit.csv --json audit.json
```

[`examples/alerts.csv`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/alerts.csv) is an 18-pattern catalogue with planted defects — an unparseable SMARTS, an exact duplicate, an aromatic/Kekulé pair, a never-matching `[#118]`, a `c[H]` that is dead without explicit hydrogens, an over-broad `[#6]`, a recursive pattern and several patterns provably contained in others ([Example data](../getting-started/example-data.md)). Against the bundled reference set it reports 1 unparseable, 1 needing explicit H, 4 over-broad, 4 that never fire (1 with a never-matching atom), 1 exact duplicate, 2 library-equivalent, 9 strictly subsumed and, with `--checks all,proof`, 12 proven redundant.

About a minute for a thousand patterns on a laptop. Two tables come back.

### Table 1: defects per pattern

```
Check                    Patterns   Share   Meaning
Unparseable                     0    0.0%   RDKit cannot compile the SMARTS
Need explicit H               358   28.6%   Contain a hydrogen atom; silently dead under the default preparation
Over-broad (>10%)               4    0.3%   Match a large share of the reference set
Never fire                    728   58.2%   Zero hits on the reference set
  ↳ never-matching atom         22    1.8%   A query atom matches no real atom — likely defective
Exact duplicates              112    9.0%   Identical SMARTS string appears earlier
Library-equivalent            210   16.8%   Identical hit set to another pattern
Strictly subsumed             354   28.3%   Another pattern's hits contain all of these
Preparation-sensitive         257   20.5%   Hit count changes with molecule preparation
No flags                       47    3.8%
```

### Table 2: reproducibility

```
Preparation   Compounds flagged   Share   Patterns firing   Verdict flips vs default
implicit-h                7,195   77.6%               523                        —
explicit-h                8,392   90.5%               519            1,527 (16.5%)
kekule                    8,207   88.5%               422            1,538 (16.6%)
```

If the flip count is a sizeable fraction of your library, the catalogue is **not reproducible unless the preparation is declared**. Pick one, record it, and pass it as `--prep` to every `filter` and `substructure` run.

## Read the per-pattern table

`audit.csv` has one row per input pattern, in input order. The `flags` column is the summary; the rest are the evidence.

| Flag | Meaning | What to do |
|---|---|---|
| `unparseable` | RDKit rejects the SMARTS | fix the syntax (`--explain` shows RDKit's message) |
| `needs-explicit-h` | contains `[H]` / `[#1]` positively | run with `--prep explicit-h`, or rewrite with H-count primitives such as `[CH2]` |
| `over-broad` | hit fraction above `--breadth-threshold` (default 10%) | is this a deliberate property filter, or a mis-specified alert? |
| `dead:never-matching-atom` | a query atom matches no atom in any reference molecule | almost certainly defective, e.g. `[N+]#[C-]` — RDKit never produces that charge pattern |
| `dead:rare-combination` | every atom is realisable; the combination does not occur in the reference set | probably fine; alert sets target rare chemistry by design |
| `dead:fires-only-with-<prep>` | matches nothing under the default preparation but does under another | declare that preparation |
| `duplicate` | identical SMARTS string appears earlier in the file | remove, or keep deliberately for provenance |
| `equivalent` | identical hit set to another pattern *on this reference set* | inspect; may differ on other chemistry |
| `subsumed` | another pattern's hit set strictly contains this one's | redundant as a rejection rule *on this reference set* |
| `prep-sensitive` | hit count differs between preparations | declare the preparation wherever the pattern is used |

!!! danger "What the audit does **not** tell you"
    - **"Never fires" is not "broken".** The reference set is drug-like ChEMBL chemistry. An alert for a reactive group that medicinal chemists avoid will fire rarely or never, and that is the alert working. The `dead:*` verdicts triage this; `never-matching-atom` is the one that signals a real defect.
    - **Equivalence and subsumption are empirical**, computed over 9,272 molecules. Two patterns indistinguishable here could differ on other chemistry. Treat them as "indistinguishable in practice on drug-like molecules", not as a theorem.
    - **Redundancy is not always a defect.** A duplicate across published sets tells you which set flagged a compound. Report it; do not silently delete it.

## Debug one pattern

```bash
chemlitmus smartsaudit --explain 'c:1(:c:c:c(:c:c:1)-[#6]=[#7]-[#7])-[#8]-[#1]'
```

!!! tip "Quote recursive SMARTS with single quotes"
    `$(...)` is command substitution in bash and zsh. Always single-quote a SMARTS that contains `$`.

Output: the normalised SMARTS, whether it needs explicit H, whether it is recursive, hits under each preparation, example matching molecules, and for **every query atom** how many reference molecules contain an atom satisfying that primitive alone. An atom with zero is the reason the pattern can never fire.

## Audit against your own chemistry

The bundled reference set is drug-like. If your library is agrochemicals, natural products or materials, use it as the reference instead:

```bash
chemlitmus smartsaudit examples/alerts.csv --library examples/library.smi      # or your own library
```

Every empirical verdict — breadth, dead, equivalence, subsumption, sensitivity — is then relative to *your* molecules, which is usually what you want.

## Know what kind of claim each finding is

Every redundancy finding carries an **evidence label** (`audit.csv` column `evidence`):

| Label | Means |
|---|---|
| `exact text duplicate` | byte-identical SMARTS — always true |
| `identical observed hit set` | same hits *on this panel*; a different library could separate them |
| `observed hit-set containment` | hits are a subset *on this panel* |
| `static proven containment` / `static proven equivalence` | proven for every molecule, with a checkable witness |
| `not observed in this reference` | zero hits here; says nothing about chemistry outside the panel |
| `unsupported` | outside the prover's atom/bond universe (recursive SMARTS, isotopes, valence, ring size) |
| `undecided` | no witness found — **not** a refutation |

The audit also records the conditions the empirical verdicts were produced under: the **match
semantics** (preparation, hydrogen handling, chirality ignored as RDKit does by default, RDKit
version), the **catalogue metadata** you supply (`--catalogue-name`, `--catalogue-version`,
`--catalogue-source`, `--catalogue-licence`) plus a content hash and the pattern count *per rule
set*, and a **reference-panel summary** (size, median heavy atoms, ring and charge fractions,
element census) so a reader can judge whether the panel resembles their chemistry.

### Never extrapolate a zero-hit result

```bash
chemlitmus smartsaudit examples/alerts.csv --holdout vendor_library.smi --holdout natural_products.smi
```

Each holdout library is matched separately, and the report names the patterns that never fired on
the main panel but *do* fire on the holdout. "Never fires" is a statement about a panel, not about
chemistry.

## Propose a cleanup without losing provenance

```bash
chemlitmus smartsaudit examples/alerts.csv --checks all,proof --cleanup proposals.csv
```

Each proposal names the pattern, what covers it, the evidence class, and one of three actions:
**recommended** (backed by a containment *proof* or a byte-identical duplicate — removal cannot
change a verdict), **review** (covered on this panel only; a different library might separate
them), or **keep provenance** (the only cover lies in a *different published set*, so removing the
rule would destroy the record of which catalogue flagged a compound). Nothing is ever removed
automatically; `--allow-cross-set-cleanup` opts into the provenance-losing case deliberately.

## Prove redundancy instead of measuring it

`smartsaudit`'s `subsumed_by` is empirical. Add the static prover to the audit, or run it alone:

```bash
chemlitmus smartsaudit examples/alerts.csv --checks all,proof -o audit.csv   # adds proven_subsumed_by / proof_status
chemlitmus smartsproof examples/alerts.csv -o proofs.csv                      # proofs only, no library needed
chemlitmus smartsproof --subsumes '[Cl][CX4]' '[#6][Cl]'             # one pair, with the witness
chemlitmus smartsproof --satisfiable '[R0;x2]'                       # can this ever match?
```

A proof holds for every molecule and comes with a checkable witness; *no witness* is undecided,
never a refutation. See [Static SMARTS containment proofs](../concepts/smarts-proofs.md) for
the method, its limits, and the result on the ChEMBL alert collection (395 of 1,251 alerts
provably redundant, 216 not analysable).

## Compare two versions of a catalogue

A catalogue changes — a new release, a vendor's re-implementation, your own clean-up after an
audit. `smartsdiff` tells you what changed in *behaviour*, not just in text:

```bash
chemlitmus smartsdiff examples/alerts.csv examples/alerts_v2.csv -o diff.csv
chemlitmus smartsdiff pains_table.csv rdkit:PAINS          # a file against RDKit's built-in catalogue
```

[`examples/alerts_v2.csv`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/examples/alerts_v2.csv) is a clean-up of `alerts.csv`: six patterns are gone (the over-broad `[#6]` and `c1ccccc1`, the duplicate, the Kekulé twin, the broken pattern and the never-matching `[#118]`), `aldehyde` and `michael_acceptor` were broadened, `phenol` and `thiol` were rewritten (`thiol` now requires a carbon neighbour, which changes nothing on this panel), and `epoxide` is new. The diff pairs 12 patterns by name and reports 4 texts rewritten, 10 with identical hits, 6 removed, 2 broadened and 1 added — and that 2,494 of the 9,272 reference molecules (26.9%) change verdict, all of them no longer flagged.

```bash
```

Patterns are paired by name, then by identical SMARTS, then by identical hit set on the
reference library (so a renamed *and* rewritten pattern that still does the same thing is
recognised). Each pair is classified:

| Behaviour | Meaning on the reference library |
|---|---|
| `same hits` | identical hit set (the text may still differ) |
| `broadened` | B hits everything A hits, plus more |
| `narrowed` | B hits a subset of what A hits |
| `shifted` | B gains some molecules and loses others |
| `broken` / `repaired` | the pattern parses on one side only |
| `added` / `removed` | present on one side only |

The headline number is **verdict changes**: how many reference molecules are flagged by one
catalogue and not the other. That is what a screen's users experience.

### Example: two implementations of PAINS

The PAINS SMARTS distributed in the ChEMBL structural-alert table and the PAINS catalogue
built into RDKit descend from the same publication and pair 480 of 481 patterns by name.
On the bundled 9,272-molecule reference set under the default (implicit-H) preparation:

```
  Molecules flagged      348 (3.8%)  vs  442 (4.8%)
  Behaviour: 462 same hits, 17 broadened, 1 shifted, 1 removed
  Verdict changes: 168 molecules (1.8%) — 131 newly flagged by B, 37 no longer flagged
```

Seventeen table patterns that never fire under implicit hydrogens (they contain explicit `[H]`
atoms; see [Molecule preparation](../concepts/molecule-preparation.md)) fire in RDKit's version,
which was rewritten with hydrogen counts; `dyes5A(27)` matches 39 molecules in the table and 6 in
RDKit. With explicit hydrogens the difference shrinks to 79 verdict changes (0.9%). The two
"PAINS" filters are not interchangeable, and which one a paper used is rarely stated.

RDKit does not expose the SMARTS text of its catalogues, so `rdkit:` sides are compared on hits
alone (`Text: no SMARTS`). Available: `rdkit:PAINS`, `rdkit:PAINS_A/B/C`, `rdkit:BRENK`,
`rdkit:NIH`, `rdkit:ZINC`, `rdkit:CHEMBL_<set>`, `rdkit:ALL`.

## Then: screen with the preparation declared

```bash
chemlitmus filter --file examples/library.csv --rules pains --prep explicit-h --output screened.csv
chemlitmus substructure 'c[H]' --file examples/library.smi --prep explicit-h
```

The `preparation` column in the output is the record that makes the screen reproducible.

## From Python

```python
from chemlitmus import audit_smarts, explain_smarts, load_patterns, load_reference_library

patterns = load_patterns("alerts.csv")                 # [(smarts, name, rule_set), ...]
mols, source = load_reference_library()                # or load_reference_library("my_library.smi")

res = audit_smarts(patterns, library=mols, library_source=source, breadth_threshold=0.05)
print(res.n_dead, res.n_dead_never_matching_atom, res.n_subsumed, res.n_clean)
print(res.sensitivity.verdict_flips)                   # {'explicit-h': 1527, 'kekule': 1538}

defective = [p for p in res.patterns if p.dead_verdict == "never-matching atom"]
```

Full model reference: [SMARTS auditing API](../reference/python-api.md#smarts-auditing).

```python
from chemlitmus import diff_smarts
res = diff_smarts("alerts_v1.csv", "alerts_v2.csv")
res.verdict_changes, res.semantic_counts
[d for d in res.patterns if d.semantic_status == "narrowed"][0].examples_lost
```

```python
from chemlitmus import subsumes, equivalent, satisfiable, prove_catalogue
subsumes("c1ccccc1[OH]", "c[OH]").proven          # True, with .witness
equivalent("[CH3]", "[C;H3]").proven               # True
satisfiable("[R0;x2]").status                      # 'unsatisfiable'
prove_catalogue([(smarts, name), ...]).n_proven_redundant
```
