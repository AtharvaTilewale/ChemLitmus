# PubChem lookups and downloads

The only part of ChemLitmus that needs a network connection. Results are cached locally in SQLite, so a second query for the same compound is instant and offline.

## Look up one compound

```bash
chemlitmus lookup "aspirin"                       # by name
chemlitmus lookup "CC(=O)Oc1ccccc1C(=O)O"         # by SMILES
chemlitmus lookup 2244 --type cid                 # by PubChem CID
chemlitmus lookup BSYNRYMUTXBXSQ-UHFFFAOYSA-N      # by InChIKey
chemlitmus lookup "aspirin" --no-cache            # force a fresh request
```

With `--type auto` (the default), the query is routed by shape: all digits → CID; 27 characters with dashes at positions 15 and 26 → InChIKey; a valid SMILES → SMILES; anything else → name. Pass `--type` explicitly when the shape is ambiguous (a name that happens to be a valid SMILES, such as `CO`).

Returned fields: `cid`, `title`, `iupac_name`, `molecular_formula`, `molecular_weight`, `canonical_smiles`, `isomeric_smiles`, `inchi`, `inchikey`, `xlogp`, `hbond_donor_count`, `hbond_acceptor_count`.

## Batch lookup

```bash
chemlitmus batch compounds.csv --output results.xlsx --format xlsx
chemlitmus batch compounds.smi --output results.json --format json --keep-duplicates
```

Reads any supported input format, deduplicates queries unless `--keep-duplicates`, runs lookups through a thread pool with rate limiting, and writes CSV, XLSX or JSON plus a `<output>_report.log` listing every query's outcome. Progress is shown live; the cache means re-running after an interruption only fetches what is missing.

## Download structures

```bash
chemlitmus download 2244 --format sdf --3d --output-dir structures
chemlitmus download "aspirin" --format png                 # 2D depiction from PubChem
chemlitmus download --file compounds.csv --format sdf --3d --output-dir structures
```

Formats from PubChem: `sdf`, `mol`, `pdb`, `png`. Existing files are skipped (resume behaviour); `--force` overwrites.

### Offline fallback and offline-only

```bash
# Try PubChem; if a compound is not there, embed it locally from its SMILES
chemlitmus download --file compounds.csv --gen missing --3d --format sdf

# Never contact PubChem; embed every SMILES locally
chemlitmus download --file compounds.smi --gen all --3d --format sdf
```

`--gen all` is fully offline. `--gen missing` is the pragmatic default for mixed libraries where most compounds are in PubChem and a few are novel. Local generation supports `sdf`, `mol` and `pdb` (not `png`), 2D coordinates or ETKDG 3D with MMFF94/UFF optimisation.

## IUPAC names

```bash
chemlitmus iupacname "CC(=O)Oc1ccccc1C(=O)O"            # offline: InChI, InChIKey, formula, exact mass
chemlitmus iupacname "CC(=O)Oc1ccccc1C(=O)O" --online   # + PubChem preferred IUPAC name, cached
chemlitmus iupacname --file compounds.csv --online --output names.csv
```

Identifiers are always computed locally. The name requires `--online` on first request and is served from the cache thereafter.

## Rate limits and etiquette

PubChem asks for no more than five requests per second per user. ChemLitmus defaults to a 0.5 s minimum spacing (two per second) with three retries and a 10 s timeout. On a shared outbound IP you may still see HTTP 503; raise `CHEMLITMUS_RATE_LIMIT_DELAY` to 1.0 or more. See [Configuration](../getting-started/configuration.md).

## Cache management

```bash
chemlitmus status      # shows the cache path and whether it is loaded
chemlitmus init        # creates directories and the database explicitly
```

The cache is a plain SQLite file (`compound_cache` and `iupac_names` tables) at `<cache_dir>/chemlitmus.db`. Delete it to start fresh; set `CHEMLITMUS_ENABLE_CACHE=false` to bypass it entirely.

## From Python

```python
from chemlitmus import lookup, lookup_by_name, lookup_file, download_structure, get_iupac_name

c = lookup("aspirin")
print(c.cid, c.iupac_name, c.inchikey)

results = lookup_file("compounds.csv", output_format="xlsx")
download_structure(2244, format="sdf", dimension="3d", output_dir="structures")
r = get_iupac_name("CCO", use_online=True); print(r.iupac_name, r.iupac_name_source)
```
