# Bundled data

## `reference_library.smi.gz`

A reference set of 9272 drug-like small molecules used by `chemlitmus smartsaudit` as the
default population against which SMARTS patterns are evaluated (hit rates, redundancy,
dead-rule detection and preparation sensitivity). Pass `--library` to substitute your own.

| | |
|---|---|
| Molecules | 9272 (3417 approved drugs with `max_phase = 4`; 5855 small-molecule sample) |
| Source | ChEMBL_37, released 2026-05-01, retrieved via the ChEMBL web services |
| Format | gzip-compressed `.smi`: `canonical_smiles chembl_id subset`, `#` comment header |
| SMILES | RDKit canonical, re-canonicalized at build time |

### Attribution

ChEMBL data are made available by the European Bioinformatics Institute under the
[Creative Commons Attribution-ShareAlike 3.0 Unported](https://creativecommons.org/licenses/by-sa/3.0/)
licence. This file is a derived subset and is redistributed under the same licence; it is
**not** covered by the MIT licence that applies to the ChemLitmus source code.

If you use results derived from this reference set, please cite ChEMBL:

Zdrazil B, et al. *The ChEMBL Database in 2023: a drug discovery platform spanning multiple
bioactivity data types and time periods.* Nucleic Acids Res. 2024;52(D1):D1180-D1192.
