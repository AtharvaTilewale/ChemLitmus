# Configuration

ChemLitmus works with no configuration. Everything below is optional and controls where files go, how PubChem is contacted, and how much is logged.

## How settings are read

Settings are resolved in this order; the first source that defines a value wins:

1. Environment variables prefixed `CHEMLITMUS_` (case-insensitive)
2. A `.env` file in the current working directory
3. Built-in defaults

Check the effective configuration at any time:

```bash
chemlitmus status
```

## Settings

### Directories

| Setting | Environment variable | Default | Purpose |
|---|---|---|---|
| `cache_dir` | `CHEMLITMUS_CACHE_DIR` | platform user cache dir (e.g. `~/.cache/chemlitmus`) | SQLite cache of PubChem results and IUPAC names |
| `data_dir` | `CHEMLITMUS_DATA_DIR` | platform user data dir | reserved for user data |
| `log_dir` | `CHEMLITMUS_LOG_DIR` | platform user log dir | per-run `.log` files from batch commands |
| `db_name` | `CHEMLITMUS_DB_NAME` | `chemlitmus.db` | cache database filename inside `cache_dir` |

Platform directories follow the `platformdirs` convention: `~/.cache`, `~/.local/share`, `~/.local/state` on Linux; `~/Library/Caches` etc. on macOS; `%LOCALAPPDATA%` on Windows.

!!! tip "Shared or CI environments"
    Point `CHEMLITMUS_CACHE_DIR` at a project-local path (e.g. `./.chemlitmus`) so the PubChem cache travels with the project and does not leak between users on a shared machine.

### Database access

These affect `resolve`, `lookup`, `batch`, `download` (without `--gen all`) and `iupacname --online`. `pubchem_timeout` and `pubchem_retries` also bound the per-request behaviour of the ChEMBL, ChEBI and KEGG providers; `resolve --timeout` caps the whole fan-out regardless.

| Setting | Environment variable | Default | Purpose |
|---|---|---|---|
| `pubchem_base_url` | `CHEMLITMUS_PUBCHEM_BASE_URL` | `https://pubchem.ncbi.nlm.nih.gov/rest/pug` | PUG REST endpoint |
| `pubchem_timeout` | `CHEMLITMUS_PUBCHEM_TIMEOUT` | `10` | seconds per HTTP request |
| `pubchem_retries` | `CHEMLITMUS_PUBCHEM_RETRIES` | `3` | retry attempts on failure |
| `rate_limit_delay` | `CHEMLITMUS_RATE_LIMIT_DELAY` | `0.5` | minimum seconds between requests (PubChem asks for ≤5/s) |
| `max_workers` | `CHEMLITMUS_MAX_WORKERS` | CPU count | threads for batch downloads |
| `batch_size` | `CHEMLITMUS_BATCH_SIZE` | `50` | records per batch |
| `enable_cache` | `CHEMLITMUS_ENABLE_CACHE` | `true` | read/write the SQLite cache |

If you share an outbound IP with other PubChem users and see 503 responses, raise `rate_limit_delay` to `1.0` or more.

### Logging

| Setting | Environment variable | Default |
|---|---|---|
| `log_level` | `CHEMLITMUS_LOG_LEVEL` | `INFO` |
| `log_format` | `CHEMLITMUS_LOG_FORMAT` | `%(asctime)s - %(name)s - %(levelname)s - %(message)s` |

Set `CHEMLITMUS_LOG_LEVEL=DEBUG` to see every PubChem request and cache hit. Batch commands additionally write a per-run log file to `log_dir`.

RDKit's own C++ warnings are suppressed throughout. The only command that surfaces RDKit parse messages is `diagnose`, which captures and interprets them deliberately.

## Example `.env`

```ini
CHEMLITMUS_CACHE_DIR=./.chemlitmus
CHEMLITMUS_RATE_LIMIT_DELAY=1.0
CHEMLITMUS_LOG_LEVEL=WARNING
```

## Programmatic access

```python
from chemlitmus.config import settings
print(settings.cache_dir, settings.rate_limit_delay)
```

`settings` is a `pydantic-settings` model; it is read once at import and is not re-read if environment variables change afterwards.
