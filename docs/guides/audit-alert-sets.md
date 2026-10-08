# Audit a structural-alert set

**Goal:** before using a SMARTS catalogue — PAINS, Brenk, Glaxo, an in-house filter — to reject compounds, find out which of its rules are broken, dead, redundant or unstable, and under which molecule preparation it was meant to be run.

## Why this is necessary

Structural-alert catalogues are redistributed as plain text and applied with whatever molecule preparation a toolkit defaults to. Nobody validates them. When ChemLitmus audited the 1,251 ChEMBL structural alerts (PAINS, BMS, SureChEMBL, MLSMR, Dundee, Inpharmatica, LINT, Glaxo) against its 9,272-molecule reference set:

- The **same catalogue flagged 77.6%, 90.5% or 88.5%** of compounds depending only on whether molecules had implicit hydrogens, explicit hydrogens or Kekulé bonds. **16.5% of compounds changed pass/fail verdict** between implicit and explicit H.
- **358 alerts (29%) contain a hydrogen atom** and therefore match nothing under RDKit's default preparation. They are silently dead unless you add hydrogens.
- **112 are byte-identical duplicates** of another alert; **354 are empirically subsumed** by another (every molecule they flag is already flagged by a broader rule).
- Only **47 of 1,251 (3.8%)** carried no flag at all.

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
chemlitmus smartsaudit alerts.csv --output audit.csv --json audit.json
```

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
chemlitmus smartsaudit alerts.csv --library my_library.smi
```

Every empirical verdict — breadth, dead, equivalence, subsumption, sensitivity — is then relative to *your* molecules, which is usually what you want.

## Then: screen with the preparation declared

```bash
chemlitmus filter --file library.csv --rules pains --prep explicit-h --output screened.csv
chemlitmus substructure 'c[H]' --file library.smi --prep explicit-h
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
