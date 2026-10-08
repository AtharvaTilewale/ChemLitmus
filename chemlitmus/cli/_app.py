"""Typer application object, console and version callback shared by all command modules."""

from typing import Optional
import typer
from rich.console import Console
from rich.panel import Panel
from rich.text import Text

from chemlitmus import (
    __version__,
)


app = typer.Typer(
    name="chemlitmus",
    help="ChemLitmus: SMILES validation and diagnosis, standardization, molecular identity, SMARTS auditing, and multi-database lookup (PubChem, ChEMBL, ChEBI, KEGG).",
    add_completion=False,
)

console = Console()

def version_callback(value: bool) -> None:
    """Print version and exit."""
    if value:
        _print_version_panel()
        raise typer.Exit()

def _print_version_panel() -> None:
    """Render the version info panel."""
    body = Text()
    body.append("\n")
    body.append("  ChemLitmus", style="bold bright_cyan")
    body.append(f"  v{__version__}", style="bold bright_green")
    body.append("\n")
    body.append("  High-performance SMILES validation, property calculation,\n", style="dim white")
    body.append("  and PubChem lookup tool for cheminformatics.\n", style="dim white")
    body.append("\n")
    body.append("  License:     ", style="white")
    body.append("MIT", style="bold white")
    body.append("\n")
    body.append("  Author:      ", style="white")
    body.append("Atharva Tilewale", style="bold white")
    body.append("\n")
    body.append("  Repository:  ", style="white")
    body.append("https://github.com/AtharvaTilewale/ChemLitmus", style="bold bright_blue underline")
    body.append("\n")
    body.append("  PyPI:        ", style="white")
    body.append("https://pypi.org/project/ChemLitmus", style="bold cyan underline")
    body.append("\n")
    panel = Panel(
        body,
        border_style="bright_cyan",
        padding=(0, 2),
    )
    console.print(panel)

@app.callback()
def main(
    version: Optional[bool] = typer.Option(
        None,
        "--version",
        "-v",
        help="Show application version and exit.",
        callback=version_callback,
        is_eager=True,
    ),
) -> None:
    """ChemLitmus CLI entry point."""
    pass
