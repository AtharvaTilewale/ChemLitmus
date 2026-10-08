# Clean a compound library

**Goal:** take a raw collection — vendor catalogue, scraped dataset, generative-model output — and produce a file where every record parses, every compound appears once at the resolution you choose, and every decision is recorded.

**Input:** any of `.csv`, `.tsv`, `.xlsx`, `.smi`, `.sdf`, `.txt`. ChemLitmus auto-detects a column named `smiles`, `canonical_smiles`, `structure` or `compound` (case-insensitive); otherwise it uses the first column. For `.smi`, the first whitespace-separated token on each line is the SMILES.

## 1. Triage: what is broken, and why

```bash
chemlitmus diagnose --file raw.csv --output diagnosis.csv
```

The terminal summary groups failures by primary cause. Open `diagnosis.csv` and sort by `primary_category`:

| Category | Usually means | Action |
|---|---|---|
| `characters` | spaces or Unicode dashes from a spreadsheet or PDF | accept the repair |
| `parentheses`, `rings` | truncated strings (cell width limit, copy-paste) | repair is a *guess*; prefer re-exporting the source |
| `brackets` | mangled `[...]` atoms, often `[C@@H` missing `]` | repair if obvious, else drop |
| `valence` | neutral four-connected N, hypervalent C | real chemistry error upstream; drop or fix by hand |
| `aromaticity` | pyrrole `n` without `[nH]`, aromatic atom outside a ring | accept `[nH]` repair; inspect the rest |

Records with `repaired_is_valid = True` have a candidate in `repaired_smiles`. **Repairs are mechanical, not chemical** — they make a string parse, they do not know what molecule was meant. Spot-check a sample before accepting them wholesale.

!!! warning "The whitespace trap"
    `CC O` is a *valid* SMILES to RDKit: it parses as ethane and silently discards `O` as a title. `diagnose` and `validate` both flag it; a plain RDKit `MolFromSmiles` check does not. If your pipeline has ever trusted RDKit's return value alone, this has happened to you.

## 2. Standardise representation

```bash
chemlitmus standardize --file triaged.csv --output standardized.csv
```

Default steps, in order: `fragment` (keep the largest organic fragment — strips counter-ions), `neutralize` (uncharge where chemically reasonable), `tautomer` (RDKit canonical tautomer), `canonical` (canonical SMILES). Choose a subset with `--steps`:

```bash
# Keep salts but neutralise and canonicalise
chemlitmus standardize --file in.csv --steps neutralize,canonical --output out.csv
```

Use `--show-diff` on a single SMILES to see what each step changed. The output CSV has `input_smiles`, `output_smiles`, `changed` and `error` columns, so nothing disappears silently.

!!! note "Standardise or not?"
    Standardising is lossy by design: the salt form, the charge state and the tautomer are gone. If any of those matters for your task (formulation, pKa modelling, crystal forms), do **not** standardise — use `identity` instead, which keeps the record and tells you how it relates to others.

## 3. Decide what "the same compound" means

```bash
chemlitmus identity --file standardized.csv --level parent --output identity.csv
```

The first table shows how many distinct compounds you have at each level:

```
Level      Distinct compounds  Collapsed records
exact                   12,480                  0
parent                  11,902                578
tautomer                11,715                765
nostereo                11,340              1,140
skeleton                11,201              1,279
```

Pick the level that matches your question — see [Molecular identity](../concepts/molecular-identity.md) for how to choose. The groups table shows every cluster of records that are the same compound at that level and *what varies* among them (salt form, tautomer, stereochemistry). `identity.csv` carries a `group_id_<level>` column; keep one row per group id to deduplicate.

```python
import pandas as pd
df = pd.read_csv("identity.csv")
dedup = df[df.is_valid].drop_duplicates("group_id_parent")
```

## 4. Screen — with the preparation declared

```bash
chemlitmus filter --file dedup.csv --rules lipinski,veber,pains --prep implicit-h --output screened.csv
chemlitmus filter --file dedup.csv --rules pains --fail --output pains_hits.csv   # keep only the flagged ones
```

Every output row records `preparation`. Structural-alert results depend on it — see [Molecule preparation](../concepts/molecule-preparation.md) — and a screen whose preparation is not recorded cannot be reproduced.

## 5. Record provenance

Keep `diagnosis.csv`, `standardized.csv` and `identity.csv` alongside the cleaned file. Between them they answer, for any compound in the output, "where did this come from and what was done to it", which is what a reviewer — or you in six months — will ask.

## Whole pipeline as a script

```bash
#!/usr/bin/env bash
set -euo pipefail
chemlitmus diagnose     --file raw.csv            --output 01_diagnosis.csv
# (review 01_diagnosis.csv; write accepted repairs + valid records to 02_triaged.csv)
chemlitmus standardize  --file 02_triaged.csv     --output 03_standardized.csv
chemlitmus identity     --file 03_standardized.csv --level parent --output 04_identity.csv
# (deduplicate on group_id_parent -> 05_dedup.csv)
chemlitmus filter       --file 05_dedup.csv --rules lipinski,pains --prep implicit-h --output 06_screened.csv
```

The review steps are deliberately manual. ChemLitmus will tell you what is wrong and propose fixes; deciding which fixes are chemically right is still your call.
