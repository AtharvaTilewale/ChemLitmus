"""
CLI entry point for ChemLitmus.

This module is referenced in pyproject.toml as the CLI entry point.
It imports and re-exports the Typer app from cli.main.
"""

from chemlitmus.cli.main import app  # noqa: F401

__all__ = ["app"]
