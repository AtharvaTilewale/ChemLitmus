# Installation

## Requirements

- Python 3.10 or newer
- RDKit 2023.09 or newer (installed automatically as a dependency)
- Linux, macOS or Windows

No compiler is needed; every dependency ships as a binary wheel.

## Install

=== "pip"

    ```bash
    pip install chemlitmus
    ```

=== "pipx (isolated CLI)"

    ```bash
    pipx install chemlitmus
    ```

    Use this if you only want the `chemlitmus` command and do not need to import the library.

=== "uv"

    ```bash
    uv tool install chemlitmus      # CLI only
    uv add chemlitmus               # as a project dependency
    ```

=== "From source"

    ```bash
    git clone https://github.com/AtharvaTilewale/ChemLitmus.git
    cd ChemLitmus
    pip install -e ".[dev]"
    ```

    The `dev` extra adds pytest, ruff, black, isort, mypy and coverage tooling.

## Verify

```bash
chemlitmus --version
chemlitmus validate "CC(=O)Oc1ccccc1C(=O)O"
```

The second command should print a table with the canonical SMILES `CC(=O)Oc1ccccc1C(=O)O`, formula `C9H8O4` and exact mass `180.0423`. If it does, RDKit is working and you are ready.

## Optional: initialise the cache

The PubChem-backed commands (`lookup`, `batch`, `download`) cache results in a local SQLite database so repeat queries are instant and offline. Create it explicitly with:

```bash
chemlitmus init
```

This is optional — the database is created on first use — but `init` also reports where the cache, data and log directories live on your system. See [Configuration](configuration.md) to move them.

## Upgrading

```bash
pip install --upgrade chemlitmus
```

or, from within the tool:

```bash
chemlitmus update --check     # see whether a newer release exists
chemlitmus update             # upgrade (asks for confirmation)
```

## Troubleshooting

**`ImportError: No module named rdkit`** — RDKit did not install. On an unusual platform, install it first from conda-forge (`conda install -c conda-forge rdkit`) and then `pip install chemlitmus --no-deps`.

**Commands print raw `[hh:mm:ss] SMILES Parse Error` lines** — you are running an older ChemLitmus. Upgrade; RDKit's log stream is now captured and only shown by `diagnose`.

**PubChem commands fail with 503 or time out** — PubChem rate-limits aggressively. ChemLitmus already throttles to about two requests per second and retries three times; if you are behind a shared IP, raise `CHEMLITMUS_RATE_LIMIT_DELAY` (see [Configuration](configuration.md)). Everything else in the tool works without network access.
