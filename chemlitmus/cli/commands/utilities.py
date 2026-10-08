"""Status, initialisation and self-update commands."""

from pathlib import Path
import typer
from rich.table import Table
from rich.panel import Panel
from rich.text import Text

from chemlitmus import (
    __version__,
)
from chemlitmus.config import settings
from chemlitmus.cli._app import app, console


@app.command(name="init")
def init_cmd() -> None:
    """Create the cache, data and log directories and the local SQLite database.

    Running this is optional - everything is created on first use - but it reports where
    ChemLitmus will keep its files on this machine and confirms the database is writable.
    """
    from chemlitmus.core.database import DatabaseManager

    created = []
    for label, path in (("cache", settings.cache_dir), ("data", settings.data_dir), ("log", settings.log_dir)):
        existed = Path(path).exists()
        Path(path).mkdir(parents=True, exist_ok=True)
        created.append((label, str(path), "exists" if existed else "created"))
    db = DatabaseManager()
    db.init_db()
    db_path = Path(settings.cache_dir) / settings.db_name
    tbl = Table(title="[bold]ChemLitmus initialised[/bold]", show_header=True, header_style="bold magenta")
    tbl.add_column("Resource", style="cyan"); tbl.add_column("Path"); tbl.add_column("Status", justify="center")
    for label, path, status in created:
        tbl.add_row(f"{label} directory", path, status)
    tbl.add_row("database", str(db_path), "ready" if db_path.exists() else "[red]missing[/red]")
    console.print(tbl)
    console.print("[dim]Override locations with CHEMLITMUS_CACHE_DIR / CHEMLITMUS_DATA_DIR / CHEMLITMUS_LOG_DIR.[/dim]")

@app.command()
def status() -> None:
    """Display system status, cache paths, and configuration."""
    table = Table(title="ChemLitmus Status & Configuration", show_header=True, header_style="bold magenta")
    table.add_column("Property", style="cyan")
    table.add_column("Value", style="green")

    table.add_row("Version", __version__)
    table.add_row("Configuration", "Loaded")
    table.add_row("Cache Enabled", str(settings.enable_cache))
    table.add_row("Cache Dir", str(settings.cache_dir))
    table.add_row("Database Path", str(settings.db_path))
    table.add_row("Log Dir", str(settings.log_dir))
    table.add_row("Log Level", settings.log_level)
    table.add_row("PubChem Base URL", settings.pubchem_base_url)
    table.add_row("Rate Limit Delay", f"{settings.rate_limit_delay}s")
    table.add_row("Max Workers", str(settings.max_workers))
    if settings.db_path.exists():
        from chemlitmus.core.database import DatabaseManager
        stats = DatabaseManager().provider_cache_stats()
        table.add_row("Cached provider records", ", ".join(f"{k}: {v}" for k, v in sorted(stats.items())) or "none")

    console.print(table)

def _detect_install_source() -> tuple[str, str]:
    """Detect whether ChemLitmus was installed from GitHub (repo / git+url) or PyPI (pip).

    Returns:
        tuple (source_type, detail_str)
        source_type: 'git_repo' | 'git_pip' | 'pip'
    """
    import importlib.metadata
    import json
    import subprocess
    import chemlitmus

    GITHUB_REPO_URL = "https://github.com/AtharvaTilewale/ChemLitmus.git"

    # 1. Check PEP 610 direct_url.json
    try:
        dist = importlib.metadata.distribution("ChemLitmus")
        direct_url_raw = dist.read_text("direct_url.json")
        if direct_url_raw:
            info = json.loads(direct_url_raw)
            url = info.get("url", "")
            if "vcs_info" in info or "github.com" in url:
                return "git_pip", f"git+{GITHUB_REPO_URL}"
            if info.get("dir_info", {}).get("editable", False) and url.startswith("file://"):
                local_dir = Path(url.replace("file:///", "").replace("file://", ""))
                if (local_dir / ".git").exists():
                    return "git_repo", str(local_dir)
    except Exception:
        pass

    # 2. Check if running inside a Git repository work-tree
    try:
        pkg_root = Path(chemlitmus.__file__).resolve().parent
        for candidate in [pkg_root, pkg_root.parent, pkg_root.parent.parent]:
            if (candidate / ".git").exists():
                return "git_repo", str(candidate)
            try:
                res = subprocess.run(
                    ["git", "-C", str(candidate), "rev-parse", "--show-toplevel"],
                    capture_output=True,
                    text=True,
                    timeout=2,
                )
                if res.returncode == 0 and res.stdout.strip():
                    return "git_repo", res.stdout.strip()
            except Exception:
                pass
    except Exception:
        pass

    # 3. Default to PyPI / pip
    return "pip", "PyPI"

@app.command(name="update")
def update_cmd(
    yes: bool = typer.Option(False, "--yes", "-y", help="Skip confirmation prompt and upgrade immediately"),
    check: bool = typer.Option(False, "--check", "-c", help="Only check for updates, do not install"),
) -> None:
    """Check for a newer version of ChemLitmus and optionally upgrade.

    Automatically detects whether ChemLitmus was installed from GitHub or PyPI,
    and uses the appropriate upgrade method (git pull / pip install git+ / pip install).
    """
    import subprocess
    import json
    import urllib.request
    import urllib.error
    from packaging.version import Version

    source_type, source_detail = _detect_install_source()
    source_label = "GitHub (local clone)" if source_type == "git_repo" else ("GitHub (git+url)" if source_type == "git_pip" else "PyPI (pip)")

    PYPI_URL = "https://pypi.org/pypi/ChemLitmus/json"

    with console.status(f"[bold green]Checking for latest version (installed from {source_label})...[/bold green]"):
        try:
            with urllib.request.urlopen(PYPI_URL, timeout=10) as resp:
                data = json.loads(resp.read().decode())
            latest_version = data["info"]["version"]
        except urllib.error.URLError as e:
            console.print(f"[red]Network error:[/red] Could not reach PyPI. Check your connection.\n{e}")
            raise typer.Exit(code=1)
        except Exception as e:
            console.print(f"[red]Error fetching version info:[/red] {e}")
            raise typer.Exit(code=1)

    try:
        current = Version(__version__)
        latest = Version(latest_version)
    except Exception:
        current_str = __version__
        latest_str = latest_version
        up_to_date = current_str == latest_str
    else:
        up_to_date = current >= latest

    # Build info panel
    body = Text()
    body.append("\n")
    body.append("  Installed: ", style="white")
    body.append(f"v{__version__}", style="bold bright_green" if up_to_date else "bold yellow")
    body.append("\n")
    body.append("  Source:    ", style="white")
    body.append(f"{source_label}", style="bold bright_cyan")
    body.append("\n")
    body.append("  Latest:    ", style="white")
    body.append(f"v{latest_version}", style="bold bright_green")
    body.append("\n")

    if up_to_date:
        body.append("\n")
        body.append("  You are up to date!", style="bold bright_green")
        body.append("\n")
        panel = Panel(
            body,
            title="[bold bright_cyan]ChemLitmus Update Check[/bold bright_cyan]",
            border_style="bright_green",
            padding=(0, 2),
        )
        console.print(panel)
    else:
        body.append("\n")
        body.append("  A new version is available: ", style="white")
        body.append(f"v{latest_version}", style="bold bright_cyan")
        body.append("\n")
        body.append("  Repository: ", style="white")
        body.append("https://github.com/AtharvaTilewale/ChemLitmus", style="bold bright_blue underline")
        body.append("\n")
        panel = Panel(
            body,
            title="[bold bright_cyan]ChemLitmus Update Check[/bold bright_cyan]",
            border_style="yellow",
            padding=(0, 2),
        )
        console.print(panel)

        if check:
            console.print("\n[dim]Run [white]chemlitmus update[/white] to upgrade.[/dim]")
            raise typer.Exit()

        # Prompt or auto-confirm
        if not yes:
            do_upgrade = typer.confirm(
                f"\nUpgrade from v{__version__} to v{latest_version} via {source_label}?",
                default=True,
            )
        else:
            do_upgrade = True

        if do_upgrade:
            console.print(f"\n[bold green]Upgrading ChemLitmus from {source_label}...[/bold green]")

            if source_type == "git_repo":
                # Upgrade via git pull in repo directory
                repo_path = source_detail
                console.print(f"[dim]Running git pull in {repo_path}...[/dim]")
                try:
                    result = subprocess.run(
                        ["git", "-C", repo_path, "pull", "origin", "main"],
                        capture_output=True,
                        text=True,
                    )
                    if result.returncode != 0:
                        # try fallback to default pull
                        result = subprocess.run(
                            ["git", "-C", repo_path, "pull"],
                            capture_output=True,
                            text=True,
                        )
                    if result.returncode == 0:
                        console.print("[bold bright_green]Successfully pulled latest changes from GitHub![/bold bright_green]")
                        console.print(f"[dim]{result.stdout.strip()}[/dim]")
                    else:
                        console.print(f"[red]git pull failed with code {result.returncode}[/red]")
                        if result.stderr:
                            console.print(f"[dim]{result.stderr.strip()}[/dim]")
                except FileNotFoundError:
                    console.print("[red]Error:[/red] git command not found. Please pull updates manually:")
                    console.print(f"  [bold white]cd {repo_path} && git pull[/bold white]")

            elif source_type == "git_pip":
                # Upgrade via pip git URL
                git_url = "git+https://github.com/AtharvaTilewale/ChemLitmus.git"
                console.print(f"[dim]Running pip install --upgrade {git_url}...[/dim]")
                try:
                    result = subprocess.run(
                        ["pip", "install", "--upgrade", git_url],
                        capture_output=True,
                        text=True,
                    )
                    if result.returncode == 0:
                        console.print(f"[bold bright_green]Successfully upgraded from GitHub to v{latest_version}![/bold bright_green]")
                    else:
                        console.print(f"[red]pip exited with code {result.returncode}[/red]")
                        if result.stderr:
                            console.print(f"[dim]{result.stderr.strip()}[/dim]")
                except FileNotFoundError:
                    console.print("[red]Error:[/red] pip not found. Upgrade manually with:")
                    console.print(f"  [bold white]pip install --upgrade {git_url}[/bold white]")

            else:
                # Upgrade via PyPI
                console.print("[dim]Running pip install --upgrade ChemLitmus...[/dim]")
                try:
                    result = subprocess.run(
                        ["pip", "install", "--upgrade", "ChemLitmus"],
                        capture_output=True,
                        text=True,
                    )
                    if result.returncode == 0:
                        console.print(f"[bold bright_green]Successfully upgraded from PyPI to v{latest_version}![/bold bright_green]")
                    else:
                        console.print(f"[red]pip exited with code {result.returncode}[/red]")
                        if result.stderr:
                            console.print(f"[dim]{result.stderr.strip()}[/dim]")
                except FileNotFoundError:
                    console.print("[red]Error:[/red] pip not found. Upgrade manually with:")
                    console.print("  [bold white]pip install --upgrade ChemLitmus[/bold white]")
        else:
            console.print("[dim]Upgrade cancelled.[/dim]")
