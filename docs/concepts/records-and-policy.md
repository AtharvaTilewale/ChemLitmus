# Records, policy and provenance

Three ideas underpin every dataset workflow in ChemLitmus: nothing is dropped, every chemical
choice is explicit, and every result says what produced it.

## Records: nothing is dropped

`read_records()` turns any supported file into a `RecordSet` of `Record`s. Each record keeps:

* `record_id` — stable internal id (`r000000`, position-based), unique even when source ids repeat;
* `position` — the data row (or SDF molecule index), so a finding points at a line in your file;
* `source_id` — the id *you* supplied; it need not be unique and is never rewritten;
* `structure` — the original text exactly as read (or the MOL block for SDF);
* `fields` — every other column or SDF property, verbatim;
* `status` — `ok` | `empty` | `invalid` | `unsupported` | `error`, with `issues` explaining why.

```python
from chemlitmus import read_records
rs = read_records("library.csv")
rs.n_total == rs.n_ok + rs.n_empty + rs.n_invalid + rs.n_unsupported + rs.n_error   # always true
rs.reconcile()                                                                      # asserts it
```

Formats: CSV, TSV, TXT, XLSX/XLS, SMI, SDF, and Parquet (optional: `pip install 'chemlitmus[parquet]'`). Column roles (structure, id, endpoint, units,
relation, split, date, source, target) are detected from recognised headers or named explicitly.
A `.smi` line is defined to carry a name after whitespace; a structure *column* is not, so
`CC O` in a CSV is an error rather than silently becoming ethane.

!!! note "Why the old parser is still there"
    `parse_compounds_file()` returns a flat list of strings and discards ids, columns and failed
    records. It is kept for backward compatibility and is now a thin wrapper over
    `read_records()`. New code should use the record model.

## Policy: every chemical choice is a choice

`ChemicalPolicy` is validated, serialisable and hashed. It covers input representations, fragment
selection, neutralisation, tautomer/stereo/isotope handling, the identity level used for grouping
and leakage, the molecule preparation and chirality setting used for matching, the alert sets, the
fingerprint and similarity threshold, the repair policy, and severity overrides per issue code.

```python
from chemlitmus import ChemicalPolicy, resolve_policy
ChemicalPolicy.preset("conservative")   # keep everything, identity at 'exact', no repairs
ChemicalPolicy.preset("parent")         # largest organic fragment, neutralised, identity at 'parent'
resolve_policy("my_policy.json").hash   # 64-hex identity of the resolved configuration
```

Two decisions deserve emphasis:

* **Largest-fragment selection is not active-compound selection.** When a record has more than one
  substantial organic component, ChemLitmus flags it (`FRAG_MULTIPLE_ORGANIC`) for review instead
  of assuming the biggest fragment is the one you care about. The original record is preserved
  either way, and the discarded components are visible in the transformation log.
* **Repairs are candidates.** `repair.mode` is `report` by default: a mechanical repair that makes
  a string parse is reported with its exact edits, never applied silently, and never assigned a
  confidence number. The original record is always retained.

## Provenance: what produced this result

Every transformation records `before`, `after`, the operation and parameters, whether the *text*
changed and whether the *identity* changed — and, when identity changed, how (`salt, counter-ion
or charge form`, `tautomer`, `stereochemistry`), in the same language as `identity` and `diff`.

Every run writes a `manifest.json`:

```json
{
  "tool_version": "1.1.0", "command": "chemlitmus audit", "created_at": "2026-10-09T...",
  "python_version": "3.11.16", "rdkit_version": "2026.03.6", "platform": "Linux x86_64",
  "inputs": [{"path": "data.csv", "sha256": "…", "size_bytes": 1234}],
  "outputs": [...], "policy_hash": "…", "policy_name": "parent",
  "rule_catalogue": {"source": "rdkit FilterCatalog", "sets": "PAINS"},
  "output_schema_version": "1", "settings": {"identity_level": "parent", "preparation": "implicit-h"}
}
```

File names only — no private directory structure, no secrets. The same input and the same resolved
policy reproduce the same chemical decisions in the same supported environment; a different RDKit
version is visible in the manifest, so a change in results can be attributed rather than guessed at.

## Limits

* Reproducibility is **within a supported environment**. RDKit updates can change aromaticity
  perception, standardisation and descriptor values; the manifest records the version so such a
  change is detectable, but ChemLitmus does not freeze RDKit's chemistry.
* The record model preserves what the file contains. It cannot recover information the file never
  had — a missing unit, an unstated assay, an ambiguous drawing.
* `status: unsupported` means ChemLitmus declined to interpret the entry (an SMI comment line, for
  instance). It is not a claim that the chemistry is invalid.
* Policy hashes identify a configuration, not a result: the same hash with a different RDKit
  version can legitimately produce different numbers.
