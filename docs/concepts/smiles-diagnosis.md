# SMILES diagnosis

RDKit's `MolFromSmiles` returns `None` when a SMILES fails. It also writes an error to its C++ log stream, which most code suppresses, and which — even when visible — reports only the first problem the parser hit. `diagnose` replaces that with a deterministic, ordered set of checks, each pointing at a position or atom, each with a suggestion, and each followed by a repair attempt where a safe one exists.

## Design principles

**Deterministic and explainable.** Every verdict comes from a rule you can read in the source. There is no model, no training set, no probability. Given the same string, `diagnose` always says the same thing and can always say *why*.

**Located.** Every problem carries a character `position` (0-based) or an `atom_index`, and the terminal output draws a `^` marker line under the input. "Unbalanced parentheses" is a message; "`(` at position 4 is never closed" is a diagnosis.

**Ordered.** Checks run in a fixed order from lexical to chemical, and a lexical finding suppresses RDKit's generic syntax error for the same string, because the lexical finding *is* the explanation.

**Repairs are mechanical, not chemical.** A repair makes a string parse. It does not know what molecule was meant. Every repair is listed in `repairs_applied`, the result is re-validated, and `repaired_is_valid` says whether it worked. Treat `repaired_smiles` as a candidate to review.

## The seven checks

| Order | Category | Finds | Repair |
|---|---|---|---|
| 1 | `characters` | whitespace inside the string; non-ASCII characters (en dashes, non-breaking spaces); characters outside the SMILES alphabet | strip whitespace; map `–` `—` `−` to `-` |
| 2 | `brackets` | `[` never closed; bracket atom not matching `[isotope? symbol chirality? H-count? charge? :map?]`; unknown element symbol | — |
| 3 | `parentheses` | `)` with no open branch; `(` never closed, each with position | remove stray `)`; remove a dangling `(` with nothing after it, else append `)` |
| 4 | `rings` | ring-closure digit opened and never closed; digit appearing before any atom; digit closing onto its own atom | remove the unclosed digit |
| 5 | `syntax` | anything RDKit's parser still rejects, with RDKit's own reported position | — |
| 6 | `valence` | explicit valence exceeding the element's permitted values, with atom index | — (suggests `[N+]` for 4-connected N, etc.) |
| 7 | `aromaticity` | ring with no Kekulé form (classic: pyrrole `n` missing `[nH]`); lower-case atom outside any ring | write one `n` as `[nH]`; upper-case a non-ring aromatic atom |

Checks 1–4 are a tokeniser and bracket/ring bookkeeping. Check 5 captures RDKit's parser log through a dedicated Python logger so the message and position are recovered without leaking to stderr. Checks 6–7 use `Chem.DetectChemistryProblems` on the unsanitised molecule, which reports every problem rather than stopping at the first.

## The whitespace case

```
CC O
```

This string is **valid** to RDKit. `MolFromSmiles("CC O")` returns ethane and discards `O` as a molecule title, because whitespace is the SMILES/title separator in `.smi` files. In a CSV cell that leaked a space, that is a silently corrupted record presented as a success.

`diagnose` reports it as invalid with the category `characters`, says what RDKit would have parsed, and repairs to `CCO`. `validate` and `identity` apply the same rule. This is the one case where ChemLitmus deliberately disagrees with RDKit about validity, and it is the case most likely to have bitten anyone who has loaded spreadsheet-exported SMILES.

## Exit codes and pipelines

In single mode `diagnose` exits `0` for valid, `2` for invalid — so it can guard a shell pipeline. In batch mode it exits `0` and reports per record; `--only-invalid` (default) keeps the output to the records that need attention, `--all` includes every record.

## Limits

- **Repairs guess.** Removing an unclosed ring digit or appending `)` produces *a* molecule, not necessarily *the* molecule. For truncated strings the only real fix is re-exporting the source.
- **One `[nH]` at a time.** The aromaticity repair tries writing each aromatic `n` as `[nH]` in turn and keeps the first that parses. Rings needing two hydrogens, or a charge, are reported but not repaired.
- **Element list is fixed.** Bracket-atom symbols are checked against the periodic table plus the aromatic lower-case set; exotic or placeholder symbols (`[R]`, `[X]`) are reported as unknown elements, which is usually what you want.
- **Valence suggestions are generic.** `[N+]` for a four-connected nitrogen is correct far more often than not, but it is a hint, not chemistry.
- **Position mapping for chemistry errors is atom-based.** `valence` and `aromaticity` problems are located by mapping RDKit's atom index back to the token that produced it; for strings that RDKit re-ordered or that contain unusual tokens the position may point at the start of the atom rather than the exact character.
