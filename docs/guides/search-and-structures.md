# Search, fingerprints and 3D structures

The everyday RDKit operations, packaged with consistent input handling, CSV output and typed return values. All offline.

## Similarity search

```bash
chemlitmus similar "CC(=O)Oc1ccccc1C(=O)O" --file library.csv --threshold 0.5 --top 10
chemlitmus similar "CCO" --file library.smi --fp-type ecfp6 --bits 4096 --output hits.csv
```

Tanimoto similarity on the chosen fingerprint, ranked, with threshold and top-N. Output columns: `rank`, `query`, `hit`, `similarity`, `fingerprint_type`.

| `--fp-type` | Description | Bits |
|---|---|---|
| `ecfp4` (default) | Morgan radius 2, connectivity | `--bits` (2048) |
| `ecfp6` | Morgan radius 3 | `--bits` |
| `fcfp4` | Morgan radius 2, pharmacophoric features | `--bits` |
| `maccs` | MACCS 166 keys | fixed 167 |
| `rdkit` | RDKit topological (Daylight-like) | `--bits` |
| `atompair` | atom pairs | `--bits` |
| `torsion` | topological torsions | `--bits` |

## Fingerprints as data

```bash
chemlitmus fingerprint "CCO" --type ecfp4
chemlitmus fingerprint --file library.csv --type all --output fp.csv
```

Single mode prints the bit count, density and hex string. Batch mode writes one row per compound with a hex-encoded fingerprint column per type. Decode in Python with `int(hexstr, 16)` or `bytes.fromhex`.

## Substructure search

```bash
chemlitmus substructure "[NX3;H2]" --file library.csv --output amines.csv
chemlitmus substructure "c1ccccc1" --file library.csv --smiles-query
chemlitmus substructure 'c[H]' --file library.smi --prep explicit-h
```

The query is SMARTS by default (`--smiles-query` to treat it as an exact SMILES fragment). Output: matching `smiles`, `match_indices` (atom indices of the first match), and `preparation`.

!!! note "`--prep` matters for substructure too"
    A pattern containing `[H]` matches nothing unless the library molecules have explicit hydrogens; a Kekulé pattern like `C1=CC=CC=C1` matches nothing against aromatic-perceived molecules. `--prep explicit-h` / `--prep kekule` prepare the library accordingly, and the choice is recorded in every output row. See [Molecule preparation](../concepts/molecule-preparation.md).

## Scaffolds and R-groups

```bash
chemlitmus scaffold "CC(=O)Oc1ccccc1C(=O)O"            # Murcko scaffold: c1ccccc1
chemlitmus scaffold --file library.csv --output scaffolds.csv
```

Acyclic molecules have no Murcko scaffold; they are reported as `acyclic` rather than as an empty string.

```bash
chemlitmus rgroup "c1ccccc1[*:1]" --file analogues.smi
chemlitmus rgroup --core "c1ccc([*:1])cc1[*:2]" --smiles "Cc1ccccc1,CCc1ccccc1O"
```

R-group decomposition against a core SMARTS with labelled attachment points `[*:1]`, `[*:2]`, …. Molecules that do not contain the core are reported, not dropped.

## Tautomers and stereochemistry

```bash
chemlitmus tautomers "Oc1nc(O)c2nc[nH]c2n1" --max 200
chemlitmus tautomers --file ligands.csv --output tautomers.csv     # one row per tautomer
```

RDKit's `TautomerEnumerator`; the canonical tautomer is marked. Use before docking when protonation/tautomer state matters.

```bash
chemlitmus stereo "C[C@H](N)C(=O)O"
chemlitmus stereo "CC(N)C(=O)O" --chiral-flag      # exit 1: a stereocentre is unassigned
```

Lists each stereocentre with its CIP label, or `?` when unassigned. `--chiral-flag` turns it into a CI gate.

## 2D and 3D structures

```bash
# From SMILES, offline, with MMFF94 (UFF fallback) optimisation
chemlitmus download "CC(=O)Oc1ccccc1C(=O)O" --gen all --3d --format sdf --output-dir structures
chemlitmus download "c1ccccc1" --gen all --2d --format mol

# Whole file, offline
chemlitmus download --file library.smi --gen all --3d --format pdb --output-dir models

# Conformer ensembles (ETKDG v3 + MMFF94), multi-model SDF
chemlitmus conformers "CC(=O)Oc1ccccc1C(=O)O" --num 50 --output aspirin_confs.sdf
```

`--gen all` never contacts PubChem — a valid SMILES is embedded directly. Formats: `sdf`, `mol`, `pdb`. For PubChem-sourced structures (including PNG depictions), omit `--gen` or use `--gen missing` to fall back to local generation only when PubChem has no record; see [PubChem lookups](pubchem.md).

## Reactions, atom maps, augmentation

```bash
chemlitmus reaction "CC(=O)O.OCC>>CC(=O)OCC"        # reactant / agent / product counts
chemlitmus atommap "CCO"                              # [CH3:1][CH2:2][OH:3]
chemlitmus augment "CCO" --num 10                     # randomised non-canonical SMILES for ML augmentation
```

## From Python

```python
from chemlitmus import (compute_similarity, compute_fingerprint, substructure_search,
                        extract_scaffold, rgroup_decomposition, enumerate_tautomers,
                        analyze_stereochemistry, generate_structure, generate_conformers)

hits = compute_similarity("CCO", library, fp_type="ecfp4", threshold=0.4, top_n=20)
fp   = compute_fingerprint("CCO", fp_type="maccs")
subs = substructure_search("[NX3;H2]", library, preparation="implicit-h")
generate_structure("CCO", "ethanol.sdf", format="sdf", dimension="3d")
```

Reference: [Python API](../reference/python-api.md).
