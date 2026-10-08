# Data and licensing

## Source code

ChemLitmus is released under the [MIT License](https://github.com/AtharvaTilewale/ChemLitmus/blob/main/LICENSE). You may use, modify and redistribute it, including commercially, provided the copyright notice is retained.

## Bundled reference library

`chemlitmus/data/reference_library.smi.gz` is the default reference population for `smartsaudit`. It is **not** covered by the MIT licence.

| | |
|---|---|
| Molecules | 9,272 — 3,417 approved drugs (ChEMBL `max_phase = 4`) plus 5,855 sampled small molecules |
| Source | ChEMBL release 37 (1 May 2026), retrieved via the ChEMBL web services |
| Format | gzip-compressed `.smi`: `canonical_smiles chembl_id subset`, with a `#` comment header |
| Processing | SMILES re-canonicalised with RDKit at build time; no other modification |
| Licence | [Creative Commons Attribution-ShareAlike 3.0 Unported](https://creativecommons.org/licenses/by-sa/3.0/) |

ChEMBL data are provided by the European Bioinformatics Institute under CC BY-SA 3.0. The bundled file is a derived subset and is redistributed under the same licence. If you redistribute ChemLitmus, or publish results derived from the reference set, you must retain the attribution and licence notice, and derivative data must be shared alike.

### Citing ChEMBL

> Zdrazil B, Felix E, Hunter F, et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. *Nucleic Acids Research*. 2024;52(D1):D1180–D1192. doi:10.1093/nar/gkad1004

### Using your own reference set

Pass `--library PATH` to `smartsaudit` (or `library=` to `audit_smarts`) to replace the bundled set entirely. No ChEMBL data are then involved, and the CC BY-SA terms do not apply to your results.

## PubChem

`lookup`, `batch`, `download` and `iupacname --online` query the [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) service. PubChem data are in the public domain in the United States; PubChem asks users to respect its [usage policy](https://pubchem.ncbi.nlm.nih.gov/docs/programmatic-access) — no more than five requests per second, and no more than 400 requests per minute. ChemLitmus throttles to two per second by default. Cached results in your local SQLite database are your own copy.

## Third-party components

ChemLitmus depends on RDKit (BSD-3-Clause), pandas (BSD-3-Clause), NumPy (BSD-3-Clause), Typer and Rich (MIT), Pydantic (MIT), requests and aiohttp (Apache-2.0), tqdm (MIT/MPL-2.0), platformdirs (MIT) and packaging (Apache-2.0/BSD). Their licences permit redistribution alongside MIT code.

## No warranty

ChemLitmus is research software. Its outputs — validity verdicts, identity groupings, alert audits, repaired SMILES — are computed heuristically and are provided without warranty of any kind. They are inputs to a scientist's judgement, not substitutes for it.
