# Contributing

Bug reports, feature requests and pull requests are welcome at [github.com/AtharvaTilewale/ChemLitmus](https://github.com/AtharvaTilewale/ChemLitmus).

## Reporting a bug

Open an issue with: the ChemLitmus and RDKit versions (`chemlitmus --version`, `python -c "import rdkit; print(rdkit.__version__)"`), the exact command or call, the input that triggers it (a single SMILES or a minimal file), what you expected and what happened. For a SMILES that behaves unexpectedly, the output of `chemlitmus diagnose` on it is usually the fastest diagnostic.

## Development setup

```bash
git clone https://github.com/AtharvaTilewale/ChemLitmus.git
cd ChemLitmus
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Tests

```bash
pytest -m "not integration"          # offline suite, ~10 minutes
pytest tests/test_smartsaudit.py     # one module, ~1 s
pytest -m integration                # PubChem-backed tests; needs network
```

Every bug fix should come with a regression test; every new command with CLI tests via `typer.testing.CliRunner` (see `tests/test_cli.py`). The offline suite must pass on Python 3.10, 3.11 and 3.12 — CI runs all three.

## Code style

```bash
ruff check chemlitmus
black chemlitmus tests
isort chemlitmus tests
```

Line length is 100. Public functions return Pydantic models and report failures as values (`error`, `is_valid=False`) rather than raising; only programming errors (unknown option names) raise `ValueError`.

## Adding a command

1. Core logic in `chemlitmus/core/<name>.py`: a function returning a Pydantic model, importing RDKit defensively behind an `_RDKIT_AVAILABLE` flag.
2. Export from `chemlitmus/core/__init__.py` and `chemlitmus/__init__.py` (add to `__all__`).
3. CLI in `chemlitmus/cli/main.py`: single-SMILES and `--file` modes, `--output` CSV, `escape()` on any SMILES/SMARTS printed through Rich.
4. Tests, then documentation: `docs/reference/cli.md`, `docs/reference/python-api.md`, and a guide or concept page if the feature needs explaining.
5. A line in `CHANGELOG.md`.

## Documentation

```bash
pip install mkdocs mkdocs-material
mkdocs serve        # live preview at http://127.0.0.1:8000
```

Documentation claims about behaviour should be verified against the tool, not inferred from the code. If a page says a command prints `valid=6 invalid=3`, run it.

## Releasing

Bump `version` in `pyproject.toml`, `__version__` in `chemlitmus/__init__.py`, and `version` in `CITATION.cff`; update `CHANGELOG.md`; tag `vX.Y.Z`; publish a GitHub release. The `publish.yml` workflow builds and uploads to PyPI via trusted publishing.
