# Contributing

Thank you for your interest in contributing to **ChemLitmus**!

We welcome bug reports, feature requests, documentation improvements, and code contributions.

---

## Reporting Issues

If you encounter a bug or have a feature request, please open an issue on GitHub and include:

- ChemLitmus version
- Python version
- Operating system
- Steps to reproduce
- Error messages (if any)

---

## Development Setup

Clone the repository:

```bash
git clone https://github.com/AtharvaTilewale/ChemLitmus.git
cd ChemLitmus
```

Create a virtual environment:

```bash
python -m venv .venv
```

Activate it:

### Windows

```bash
.venv\Scripts\activate
```

### Linux/macOS

```bash
source .venv/bin/activate
```

Install in editable mode:

```bash
pip install -e ".[dev]"
```

---

## Project Structure

```
ChemLitmus/
├── .github/workflows/      # CI/CD workflows
├── chemlitmus/          # Main package
│   ├── __init__.py         # Public API
│   ├── config.py           # Configuration management
│   ├── logging_config.py   # Logging setup
│   ├── cli.py              # CLI entry point
│   ├── cli/
│   │   ├── __init__.py
│   │   └── main.py         # Typer CLI application
│   ├── core/               # Core functionality
│   │   ├── __init__.py
│   │   ├── smiles.py       # SMILES validation
│   │   ├── pubchem.py      # PubChem API
│   │   └── database.py     # SQLite caching
│   └── utils/              # Utilities
│       ├── __init__.py
│       ├── file_io.py      # File parsing
│       ├── export.py       # Export formats
│       └── parsers.py      # Input parsers
├── tests/                  # Test suite
├── docs/                   # Documentation
├── pyproject.toml          # Package metadata & dependencies
├── README.md               # This file
└── LICENSE                 # MIT License
```

## Running Tests

Execute the test suite with:

```bash
pytest
```

---

## Coding Guidelines

- Follow PEP 8 style guidelines.
- Add docstrings for public functions.
- Keep functions modular and well documented.
- Include tests for new features whenever possible.

---

## Pull Requests

Before submitting a pull request:

- Create a feature branch.
- Keep changes focused on a single feature or fix.
- Update documentation if necessary.
- Ensure all tests pass.

Example workflow:

```bash
git checkout -b feature/my-feature
git commit -m "Add my feature"
git push origin feature/my-feature
```

Then open a Pull Request on GitHub.

---

## Questions

For questions or discussions, please use the GitHub Discussions page.

Thank you for helping improve **ChemLitmus**!