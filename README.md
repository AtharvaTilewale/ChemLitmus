# ChemLitmus

<p align="center">
  <img src="https://raw.githubusercontent.com/AtharvaTilewale/ChemLitmus/main/docs/assets/chemlitmus-logo.png" alt="ChemLitmus" width="620" />
</p>

A high-performance, production-grade tool for SMILES validation, PubChem lookup, and chemical structure retrieval.

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Features

- **SMILES Validation & Canonicalization** - Validate and standardize SMILES strings using RDKit
- **Multi-format Input** - Support for CSV, TSV, XLSX, SMI, SDF, and TXT files
- **Smart Auto-detection** - Automatically identify SMILES columns
- **PubChem Lookup** - Search by SMILES, CID, Name, InChI, and InChIKey
- **Rich Metadata** - Retrieve IUPAC name, molecular formula, mass, descriptors
- **Structure Downloads** - Get 2D/3D SDF, MOL, PDB, and PNG formats from PubChem
- **Offline Molecule Generation (--gen)** - Generate 2D and 3D conformations (SDF, MOL, PDB) offline from SMILES via RDKit with forcefield optimization
- **Molecular Fingerprints** - ECFP4, ECFP6, FCFP4, MACCS, RDKit, AtomPair, Torsion fingerprints offline (`fingerprint` command)
- **Similarity Search** - Tanimoto-based library search with threshold and top-N ranking (`similar` command)
- **Drug-Likeness Filtering** - Lipinski Ro5, Veber, Ghose, Egan, Ro3, PAINS alerts, and QED scoring (`filter` command)
- **SMARTS Pattern Auditing** - Validate structural-alert and substructure-filter sets for unparseable, dead, redundant, over-broad and preparation-sensitive patterns (`smartsaudit` command)
- **Layered Molecular Identity** - Group a collection at exact / parent / tautomer / nostereo / skeleton / formula level and see what varies (`identity` command)
- **Structure-Aware Library Diff** - Compare two compound collections by chemical identity, not SMILES text: added, removed, changed-with-reason (`diff` command)
- **SMILES Diagnosis** - Explain *where* and *why* a SMILES fails, with character positions, suggested fixes and safe mechanical repairs (`diagnose` command)
- **Batch Processing** - Process hundreds of compounds with progress tracking
- **Async/Multithreading** - Fast parallel downloads with retry logic
- **Caching** - SQLite database for storing results locally
- **Multiple Exports** - Save results as CSV, Excel, or JSON
- **Python API** - Use directly in your scripts via chemlitmus module
- **CLI Tool** - Full-featured command-line interface with chemlitmus command

## Installation

### From PyPI 

```bash
pip install chemlitmus
```

### Development Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/AtharvaTilewale/ChemLitmus.git
cd ChemLitmus
pip install -e ".[dev]"
```

## Quick Start

### CLI Usage

```bash
# Show configuration and status
chemlitmus status

# Initialize directories and database
chemlitmus init

# Lookup a single compound (by SMILES, CID, or Chemical Name)
chemlitmus lookup "c1ccccc1"  # Benzene
chemlitmus lookup "aspirin"
chemlitmus lookup 5282253 --type cid

# Batch process a file to retrieve metadata
chemlitmus batch compounds.csv --output results.xlsx --format xlsx

# Download structure from PubChem
chemlitmus download 5282253 --format sdf --3d

# Generate 3D structure offline from SMILES using RDKit (--gen all)
chemlitmus download "CC(=O)OC1=CC=CC=C1C(=O)O" --gen all --3d --format sdf

# Generate 2D MOL structure locally from SMILES
chemlitmus download "c1ccccc1" --gen all --2d --format mol

# Batch download with offline fallback for missing structures (--gen missing)
chemlitmus download --file compounds.csv --gen missing --3d --format sdf --output-dir ./structures/

# Batch generate all structures offline from a SMILES file (--gen all)
chemlitmus download --file compounds.smi --gen all --3d --format pdb --output-dir ./3d_models/
```

#
### Molecular Fingerprints

```bash
# Generate ECFP4 fingerprint (single compound)
chemlitmus fingerprint "CC(=O)OC1=CC=CC=C1C(=O)O" --type ecfp4

# All 7 fingerprint types at once
chemlitmus fingerprint "CC(=O)OC1=CC=CC=C1C(=O)O" --type all

# Batch - save to CSV
chemlitmus fingerprint --file compounds.smi --type maccs --output fingerprints.csv
```

### Similarity Search

```bash
# Top-10 most similar compounds (Tanimoto >= 0.5)
chemlitmus similar "CC(=O)OC1=CC=CC=C1C(=O)O" --file library.smi --threshold 0.5 --top 10

# Save hits to CSV
chemlitmus similar "CCO" --file compounds.csv --fp-type ecfp6 --output hits.csv
```

### Drug-Likeness Filtering

```bash
# Single compound - all rules
chemlitmus filter "CC(=O)OC1=CC=CC=C1C(=O)O"

# Batch - keep only Lipinski-compliant, PAINS-free compounds
chemlitmus filter --file compounds.csv --rules lipinski,veber,pains --output drug_like.csv

# Remove PAINS compounds from a library
chemlitmus filter --file library.csv --rules pains --output no_pains.csv

# Lead-like compounds with minimum QED 0.5
chemlitmus filter --file library.csv --rules ro3 --qed-min 0.5 --output leads.csv
```

## Python API

```python
from chemlitmus import lookup, lookup_file, download_structure, generate_structure, validate_smiles

# Lookup single compound
result = lookup("c1ccccc1")
print(result.cid, result.iupac_name)

# Process batch file
results = lookup_file("compounds.csv", output_format="xlsx")

# Download structure from PubChem
download_structure(5282253, format="sdf", dimension="3d")

# Generate 2D or 3D structure offline from SMILES
generate_structure(
    smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
    output_path="aspirin_3d.sdf",
    format="sdf",
    dimension="3d",
    title="Aspirin"
)
``` 
For more detailed API documentation, see the **[API Reference](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/docs/api_reference.md)** page.

## Documentation

For complete tutorials and advanced usage examples, see the **[Practical Guide](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/docs/practical_guide.md)** or visit the **[official documentation](https://chemlitmus.readthedocs.io/)** on Read the Docs.

## Requirements

- Python 3.10+
- RDKit (cheminformatics library)
- pandas (data handling)
- requests/aiohttp (HTTP)
- typer (CLI framework)
- rich/tqdm (UI/progress)

## Configuration

For configuration and architecture details, see the **[Configuration & Architecture](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/docs/configuration.md)** page.


## Standardization Pipeline

Clean up "dirty" SMILES: strip salts, neutralize charges, canonicalize tautomers.
Fully **offline** — powered by RDKit `MolStandardize`.

```bash
# Strip salt, neutralize, and canonicalize in one go
chemlitmus standardize "[Na+].[OH-].CC(=O)[O-]"
# Output: CC(=O)O

# Show step-by-step diff
chemlitmus standardize "[Na+].[OH-].CC(=O)[O-]" --show-diff

# Custom steps
chemlitmus standardize "O=C([O-])c1ccccc1" --steps neutralize,canonical

# Batch from CSV
chemlitmus standardize --file dirty.csv --output clean.csv
```

### Python API

```python
from chemlitmus import standardize_smiles

result = standardize_smiles("[Na+].[OH-].CC(=O)[O-]")
print(result.output_smiles)    # CC(=O)O
print(result.changed)          # True

# Inspect what changed per step
for s in result.step_results:
    if s.changed:
        print(f"{s.step}: {s.input_smiles} -> {s.output_smiles}")
```

---


## Tautomer Enumeration

Enumerate all plausible tautomers for a given SMILES string. Critical for protein-ligand docking preparation where different tautomeric states have different binding affinities.

```bash
# Single compound
chemlitmus tautomers "Oc1nc(O)c2nc[nH]c2n1"

# Batch from CSV (outputs one row per tautomer)
chemlitmus tautomers --file ligands.csv --output tautomers.csv
```

### Python API

```python
from chemlitmus import enumerate_tautomers

result = enumerate_tautomers("Oc1nc(O)c2nc[nH]c2n1", max_tautomers=1000)
print(result.num_tautomers)        # e.g., 15
print(result.canonical_tautomer)   # O=c1[nH]c(=O)c2[nH]cnc2[nH]1

for t in result.tautomers:
    print(t)
```

---

## IUPAC Identifier Generation

Generate InChI, InChIKey, formula, MW offline. Fetch IUPAC systematic name via PubChem with local caching.

```bash
# Offline: InChI, InChIKey, formula, MW
chemlitmus iupacname "CC(=O)OC1=CC=CC=C1C(=O)O"

# With IUPAC name from PubChem (cached after first fetch)
chemlitmus iupacname "CCO" --online

# Batch
chemlitmus iupacname --file compounds.smi --online --output identifiers.csv
```

### Python API

```python
from chemlitmus import get_iupac_name

result = get_iupac_name("CCO", use_online=True)
print(result.inchi)              # InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3
print(result.inchikey)           # LFQSCWFLJHTTHZ-UHFFFAOYSA-N
print(result.molecular_formula)  # C2H6O
print(result.iupac_name)         # ethanol
print(result.iupac_name_source)  # pubchem (or 'cache' on repeat calls)
```

## Molecular Identity and Library Comparison

The same compound is written many ways: as a salt or a free base, as either tautomer, with or
without stereo. `identity` computes nested identity keys and groups a collection at the level you
care about; `diff` uses the same keys to compare two collections the way a chemist would.

| Level | Identical when |
|---|---|
| `exact` | canonical SMILES match, salts and charges included |
| `parent` | largest fragment, neutralised, stereo retained |
| `tautomer` | parent, tautomer-insensitive (stereo retained) |
| `nostereo` | parent, stereochemistry removed |
| `skeleton` | parent, no stereo, tautomer-insensitive |
| `formula` | same molecular formula of the parent |

```bash
# Keys for one molecule at every level
chemlitmus identity "CC(N)C(=O)O.[Na+].[Cl-]"

# How many distinct compounds does my file really contain, and what varies among duplicates?
chemlitmus identity --file library.csv --level parent --output identity.csv

# Compare two releases / vendor catalogues / generated sets
chemlitmus diff release_34.smi release_35.smi --level parent --output diff.csv
chemlitmus diff old.csv new.csv --level skeleton --json diff.json --include-unchanged
```

`diff` reports compounds **added**, **removed**, **unchanged**, and **changed** — present in both but
written differently — with the reason: *salt, counter-ion or charge form*, *tautomer*,
*stereochemistry*, or combinations. On the bundled reference set, 3,417 approved-drug records
collapse to 2,289 distinct parent compounds; the other 1,128 are salt forms of another record.

```python
from chemlitmus import compute_identity, group_by_identity, diff_libraries, describe_difference

k = compute_identity("C[C@H](N)C(=O)O")
print(k.parent, k.nostereo, k.skeleton)

rep = group_by_identity(smiles_list, level="skeleton")
for g in rep.groups:                       # multi-member groups, largest first
    print(g.size, g.differs_by, g.smiles)

d = diff_libraries(old_smiles, new_smiles, level="parent")
print(d.n_added, d.n_removed, d.n_changed, d.changes_by_kind)
print(describe_difference(compute_identity("Oc1ccccn1"), compute_identity("O=c1cccc[nH]1")))  # tautomer
```

---

## SMILES Diagnosis

RDKit says a SMILES failed. `diagnose` says where and why, in a fixed order of checks, each
pointing at a character position or an atom, and attempts the mechanical repairs that are safe:

```bash
chemlitmus diagnose "C1CC(C"
chemlitmus diagnose "c1cncc1"                 # → pyrrole nitrogen needs [nH]; repaired: c1cc[nH]c1
chemlitmus diagnose "CC O"                    # → RDKit silently parsed only 'CC'; whitespace splits the record
chemlitmus diagnose --file generated.smi --output diagnosis.csv
```

| Check | Finds | Repair attempted |
|---|---|---|
| characters | whitespace (which RDKit silently treats as end-of-SMILES), non-ASCII dashes, invalid symbols | strip / normalise |
| brackets | malformed `[...]` atoms, unknown element symbols | — |
| parentheses | unmatched `(` or `)`, with position | close or remove dangling branch |
| rings | ring digits never closed, or closed onto the same atom | drop the unclosed digit |
| syntax | anything RDKit's parser still rejects, with its reported position | — |
| valence | atoms over their permitted valence, with atom index | — (suggests `[N+]` etc.) |
| aromaticity | rings with no Kekulé form; aromatic atoms outside rings | `[nH]`; upper-case the stray atom |

Repairs are mechanical, not chemical: every change is listed in `repairs_applied` and the result is
re-validated. Exit code is 2 for an invalid single SMILES, so it works as a guard in shell pipelines.

```python
from chemlitmus import diagnose_smiles
d = diagnose_smiles("c1cc2ccccc2n1c")
print(d.is_valid, [(p.category, p.position, p.message) for p in d.problems])
print(d.repaired_smiles, d.repaired_is_valid, d.repairs_applied)
```

---

## SMARTS Pattern Auditing

Structural-alert sets (PAINS, Brenk, Glaxo, in-house filters) are passed around as SMARTS text and applied
with whatever molecule preparation a toolkit happens to use. Nobody checks them. `smartsaudit` evaluates a
pattern set against a reference population of real molecules and reports what is actually wrong with it.
Fully **offline**; matching runs multithreaded in C++ through RDKit's `SubstructLibrary`.

```bash
# Audit a catalogue (CSV with a 'smarts' column, or one SMARTS per line)
chemlitmus smartsaudit alerts.csv

# Write the per-pattern table and the full result
chemlitmus smartsaudit alerts.csv --output audit.csv --json audit.json

# Only the checks you want; tighten the breadth threshold to 5 %
chemlitmus smartsaudit alerts.csv --checks redundancy,dead --breadth-threshold 0.05

# Audit against your own compound collection instead of the bundled reference set
chemlitmus smartsaudit alerts.csv --library my_library.smi

# Debug one pattern: atom-by-atom realisability, hits per preparation, example matches
chemlitmus smartsaudit --explain "c:1(:c:c:c(:c:c:1)-[#6]=[#7]-[#7])-[#8]-[#1]"
```

| Check | What it reports |
|---|---|
| **compile** | Patterns RDKit cannot parse; patterns that contain a hydrogen atom and are therefore silently dead unless molecules carry explicit hydrogens |
| **breadth** | Patterns matching more than a threshold share of the reference set (default 10 %) |
| **dead** | Patterns that match nothing, triaged into *never-matching atom* (a query atom no real atom satisfies — likely defective), *rare combination* (every atom is realisable), or *fires only with &lt;preparation&gt;* |
| **redundancy** | Exact duplicate SMARTS; patterns with identical hit sets; patterns strictly subsumed by another pattern in the set |
| **sensitivity** | How hit counts and per-compound pass/fail verdicts change when molecules are prepared with implicit hydrogens, explicit hydrogens, or kekulized bonds |

The bundled reference set is 9,272 ChEMBL molecules (3,417 approved drugs plus a small-molecule sample); see
[`chemlitmus/data/README.md`](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/chemlitmus/data/README.md)
for provenance and licence. A thousand-pattern catalogue audits in about a minute on a laptop.

### Declaring molecule preparation

Because alert catalogues give different answers under different preparations, `filter` and `substructure`
now take `--prep implicit-h|explicit-h|kekule` and record the choice in their output, so a screen can be
reproduced exactly:

```bash
chemlitmus filter --file library.csv --rules pains --prep explicit-h --output screened.csv
chemlitmus substructure "c[H]" --file library.smi --prep explicit-h
```

### Python API

```python
from chemlitmus import audit_smarts, explain_smarts, load_patterns, load_reference_library

patterns = load_patterns("alerts.csv")              # [(smarts, name, rule_set), ...]
mols, source = load_reference_library()             # bundled set, or pass a path

result = audit_smarts(patterns, library=mols, library_source=source)
print(result.n_dead, result.n_subsumed, result.sensitivity.verdict_flips)
for p in result.patterns:
    if not p.clean:
        print(p.name, p.flags)

e = explain_smarts("[N+]#[C-]", library=mols)
print(e.verdict)                                    # dead: never-matching atom
print([(a.query, a.n_matching_molecules) for a in e.atoms])
```

---
## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (git checkout -b feature/amazing-feature)
3. Commit changes (git commit -m 'Add amazing feature')
4. Push to branch (git push origin feature/amazing-feature)
5. Open a Pull Request

For more details, see the **[Contributing Guide](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/docs/contributing.md)**.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use ChemLitmus in your research, please cite:

```bibtex
@software{chemlitmus2026,
  author={Atharva Tilewale},
  title={ChemLitmus: A High-Performance SMILES Validation and PubChem Lookup Tool},
  version={1.0.0},
  year={2026},
  url={https://github.com/AtharvaTilewale/ChemLitmus}
}
```

## Support

- **Documentation**: [https://chemlitmus.readthedocs.io](https://chemlitmus.readthedocs.io)
- **Issues**: [https://github.com/AtharvaTilewale/ChemLitmus/issues](https://github.com/AtharvaTilewale/ChemLitmus/issues)
- **Discussions**: [https://github.com/AtharvaTilewale/ChemLitmus/discussions](https://github.com/AtharvaTilewale/ChemLitmus/discussions)

## Changelog

See [CHANGELOG.md](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/CHANGELOG.md) for version history.

---

Made with ❤️ for the cheminformatics community
