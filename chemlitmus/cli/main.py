"""ChemLitmus CLI application.

The Typer ``app`` lives in :mod:`chemlitmus.cli._app`; each module under
:mod:`chemlitmus.cli.commands` registers its commands on import. Import order
below fixes the order in which commands are listed by ``--help``.
"""

from chemlitmus.cli._app import app, console, main, version_callback  # noqa: F401
from chemlitmus.cli.commands import (  # noqa: F401
    validation,
    identity,
    smarts,
    search,
    structures,
    databases,
    utilities,
)

__all__ = ["app", "console", "main", "version_callback"]

if __name__ == "__main__":
    app()
