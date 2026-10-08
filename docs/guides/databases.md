# Look compounds up across databases

ChemLitmus talks to four public compound databases — **PubChem**, **ChEMBL**, **ChEBI** and
**KEGG** — through one interface, and uses **UniChem** to collect identifiers from more than
thirty others. You can query any one of them directly, or ask all of them at once and let
ChemLitmus reconcile the answers.

## One query, several databases

```bash
chemlitmus resolve aspirin
```

```
                       Resolve: aspirin
┏━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━┓
┃ Source  ┃ Status┃ ID          ┃ Name                 ┃ Time ┃
┡━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━┩
│ pubchem │ found │ 2244        │ Aspirin              │ 0.9s │
│ chembl  │ found │ CHEMBL25    │ ASPIRIN              │ 0.8s │
│ chebi   │ found │ CHEBI:15365 │ acetylsalicylic acid │ 1.4s │
└─────────┴───────┴─────────────┴──────────────────────┴──────┘
  sources agree on the structure  ·  InChIKey BSYNRYMUTXBXSQ-UHFFFAOYSA-N
```

Below the per-source table ChemLitmus prints a **merged record** (name, formula, masses, SMILES,
InChIKey, descriptors, synonyms — each field taken from the first source that supplies it, in the
order you listed the sources), a **cross-reference table** (CAS, DrugBank, HMDB, KEGG, PDBe,
SureChEMBL, ChemSpider, … via UniChem plus the databases' own links), and any source-specific
extras (ChEMBL's clinical phase and ATC codes, ChEBI's definition and curation stars, KEGG
pathways and enzymes, PubChem's IUPAC name).

The query can be a **name**, a **SMILES**, an **InChIKey**, or a **native identifier** — a PubChem
CID, `CHEMBL25`, `CHEBI:15365` or `C01405`. By default the type is detected; force it with
`--type name|smiles|inchikey|id`.

### Choosing sources

```bash
chemlitmus resolve caffeine -s chembl,chebi          # just these two
chemlitmus resolve caffeine -s all                   # PubChem, ChEMBL, ChEBI and KEGG
chemlitmus resolve C07481 -s kegg                    # a single database
```

The default set is `pubchem,chembl,chebi`. KEGG is opt-in because its terms of use restrict
automated access to academic, non-commercial work (see [Data and licensing](../project/data-and-licensing.md)).

### Why the agreement line matters

Each database applies its own naming rules, so the *same name* can land on *different
structures*. ChemLitmus compares the InChIKeys returned by every source:

* **agree** — all sources returned the same InChIKey.
* **disagree** — at least one source returned a different structure. The per-source table shows
  which one, so you can decide which to trust.
* **unknown** — fewer than two sources returned an InChIKey.

When one source disagrees with the majority (or does not know the name at all), ChemLitmus
re-queries it **by the consensus InChIKey**. A name that a free-text search would have matched
to a derivative (ChEBI's top hit for "aspirin", for example, is *aspirin-triggered resolvin D3*)
is corrected to the entry with the right structure, and a database that lacks the synonym is
still found through the structure.

### Timeouts and busy servers

Everything for one query runs in parallel under a single budget (`--timeout`, default 40 s).
A server that is down or busy fails fast and is reported as an **error** row; the other sources
still return. PubChem's rate limit (five requests per second) and the ChEMBL/ChEBI/KEGG service
etiquette are enforced per database.

## Batches

```bash
chemlitmus resolve --file compounds.txt -s chembl,chebi --output resolved.csv --json resolved.json
```

The input is one query per line (or a CSV/TSV with a `name`, `id`, `smiles` or `query` column).
`--output` writes one row per *(query, source)* with id, name, SMILES, InChIKey, formula,
molecular weight, URL, and the per-query agreement verdict. `--json` keeps everything,
including cross-references and source extras.

## PubChem only: `lookup` and `download`

The `lookup` and `download` commands are PubChem-specific. `lookup` returns the full PubChem
property set and caches results in the local SQLite database; `download` fetches 2D/3D structures
as SDF/MOL. See the [CLI reference](../reference/cli.md#lookup) for the options.

```bash
chemlitmus lookup aspirin
chemlitmus lookup 2244 --type cid --json
chemlitmus download CCO --format sdf --gen 3d -o ethanol.sdf
```

## From Python

```python
from chemlitmus import resolve, get_provider, unichem_xrefs

res = resolve("ibuprofen", sources=["chembl", "chebi", "pubchem"])
res.agreement                 # 'agree'
res.consensus_inchikey        # 'HEFNNWSXXWATRW-UHFFFAOYSA-N'
res.merged.formula            # 'C13H18O2'
res.cross_refs["drugbank"]    # 'DB01050'
res.by_source("chembl").extra["max_phase"]

# a single database, native schema mapped onto CompoundRecord
rec = get_provider("kegg").lookup("C01405")
rec.cross_refs                # {'cas': '50-78-2', 'pubchem_sid': '4594', 'chebi': 'CHEBI:15365', ...}

# identifiers across ~40 sources from one InChIKey
unichem_xrefs("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
```

`CompoundRecord` is the same shape for every database: `source`, `source_id`, `name`,
`synonyms`, `smiles`, `inchi`, `inchikey`, `formula`, `molecular_weight`, `monoisotopic_mass`,
`charge`, `xlogp`, `tpsa`, `hbd`, `hba`, `rotatable_bonds`, `cross_refs`, `url`, `extra`.
