"""ChemLitmus CLI Application."""

from pathlib import Path
from typing import Optional
import re
import hashlib
import typer
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text
from rich.markup import escape
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from chemlitmus import (
    __version__, lookup, download_structure, generate_structure, validate_smiles,
    compute_fingerprint, apply_filters, compute_similarity,
    substructure_search, standardize_smiles, StandardizeResult, STANDARDIZE_STEPS,
    get_iupac_name, IUPACResult,
    enumerate_tautomers, TautomerResult,
    validate_reaction, generate_conformers, extract_scaffold, analyze_stereochemistry, rgroup_decomposition, augment_smiles, map_atoms, audit_smarts, explain_smarts, load_patterns, load_reference_library,
    SmartsAuditResult, SmartsExplanation, compute_identity, group_by_identity, IDENTITY_LEVELS,
    diff_libraries, diagnose_smiles, SmilesDiagnosis,
)
from chemlitmus.config import settings
from chemlitmus.logging_config import logger
from chemlitmus.core.pubchem import PubChemCompound
from chemlitmus.utils.parsers import parse_compounds_file
from chemlitmus.utils.export import export_results

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


# =============================================================================
# resolve command (multi-database)
# =============================================================================

@app.command(name="resolve")
def resolve_cmd(
    query: str = typer.Argument(None, help="Name, SMILES, InChIKey, or native ID (CID, CHEMBL…, CHEBI:…, C…)."),
    file: Optional[Path] = typer.Option(None, "--file", "-f", help="Batch input (CSV/.smi/.txt)."),
    sources: str = typer.Option("pubchem,chembl,chebi", "--sources", "-s", help="Comma-separated: pubchem, chembl, chebi, kegg, or all."),
    query_type: str = typer.Option("auto", "--type", "-t", help="auto, name, smiles, inchikey, or id."),
    no_unichem: bool = typer.Option(False, "--no-unichem", help="Skip UniChem cross-reference lookup."),
    timeout: float = typer.Option(40.0, "--timeout", help="Seconds allowed for the whole fan-out per query."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Batch: write one row per (query, source) to CSV."),
    json_out: Optional[Path] = typer.Option(None, "--json", help="Write the complete result(s) to JSON."),
) -> None:
    """Look a compound up in several databases at once and reconcile the answers.

    Queries PubChem, ChEMBL, ChEBI and/or KEGG in parallel, returns one record per source in a
    common schema, merges them, reports whether the sources agree on the structure (by InChIKey),
    and gathers cross-database identifiers via UniChem. A name one database knows under a
    different label is recovered in the others through the structure.
    """
    import csv
    import json

    from chemlitmus.providers import PROVIDERS, resolve

    if query is None and file is None:
        console.print("[red]Error:[/red] Provide a query or --file.")
        raise typer.Exit(code=1)
    src_list = [s.strip().lower() for s in sources.split(",") if s.strip()]
    bad = [s for s in src_list if s != "all" and s not in PROVIDERS]
    if bad:
        console.print(f"[red]Error:[/red] Unknown source(s) {bad}. Available: {sorted(PROVIDERS)} or 'all'.")
        raise typer.Exit(code=1)
    if query_type.lower() not in ("auto", "name", "smiles", "inchikey", "id"):
        console.print(f"[red]Error:[/red] Unknown --type {query_type!r}. Valid: auto, name, smiles, inchikey, id.")
        raise typer.Exit(code=1)

    def _print_one(res) -> None:
        head = Table(title=f"[bold]Resolve: {escape(res.query)}[/bold]", show_header=True, header_style="bold magenta")
        head.add_column("Source", style="cyan"); head.add_column("Status"); head.add_column("ID"); head.add_column("Name", max_width=40, overflow="fold"); head.add_column("Time", justify="right")
        for o in res.outcomes:
            rec = res.by_source(o.source)
            status = {"found": "[green]found[/green]", "not found": "[yellow]not found[/yellow]"}.get(o.status, "[red]error[/red]")
            head.add_row(o.source, status, escape(rec.source_id) if rec else "", escape(rec.name or "") if rec else escape(o.error or ""), f"{o.seconds:.1f}s")
        console.print(head)
        if not res.found:
            return
        agree = {"agree": "[green]sources agree on the structure[/green]", "disagree": "[red]sources DISAGREE on the structure[/red]",
                 "unknown": "[dim]agreement unknown (fewer than two InChIKeys)[/dim]"}[res.agreement]
        console.print(f"  {agree}" + (f"  ·  InChIKey {res.consensus_inchikey}" if res.consensus_inchikey else ""))
        mg = res.merged
        t2 = Table(title="Merged record", show_header=True, header_style="bold blue")
        t2.add_column("Field", style="cyan"); t2.add_column("Value", overflow="fold")
        for label, val in (("Name", mg.name), ("Formula", mg.formula), ("MW", mg.molecular_weight), ("Monoisotopic", mg.monoisotopic_mass),
                           ("SMILES", mg.smiles), ("InChIKey", mg.inchikey), ("XLogP", mg.xlogp), ("HBD / HBA", f"{mg.hbd} / {mg.hba}" if mg.hbd is not None else None),
                           ("TPSA", mg.tpsa), ("Synonyms", ", ".join(mg.synonyms[:8]) if mg.synonyms else None)):
            if val not in (None, ""):
                t2.add_row(label, escape(str(val)))
        console.print(t2)
        if res.cross_refs:
            xr = Table(title="Cross-references", show_header=True, header_style="bold blue")
            xr.add_column("Database", style="cyan"); xr.add_column("Identifier")
            for k, v in res.cross_refs.items():
                xr.add_row(k, escape(str(v)))
            console.print(xr)
        for rec in res.records:
            if rec.extra:
                bits = []
                for k, v in rec.extra.items():
                    if v in (None, "", [], {}):
                        continue
                    sv = ", ".join(map(str, v)) if isinstance(v, list) else str(v)
                    bits.append(f"{k}={sv[:60]}")
                if bits:
                    console.print(f"  [dim]{rec.source}:[/dim] " + escape("; ".join(bits[:8])))

    if query is not None and file is None:
        with console.status(f"[bold green]Resolving '{escape(query)}' across {', '.join(src_list)}...[/bold green]"):
            res = resolve(query, sources=src_list, query_type=query_type, unichem=not no_unichem, timeout=timeout)
        _print_one(res)
        if json_out:
            json_out.parent.mkdir(parents=True, exist_ok=True)
            json_out.write_text(json.dumps(res.model_dump(), indent=2, default=str))
            console.print(f"[green]Saved:[/green] {json_out}")
        raise typer.Exit(code=0 if res.found else 1)

    queries = parse_compounds_file(file)
    results = []
    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"), BarColumn(), TaskProgressColumn(), console=console) as progress:
        task = progress.add_task("Resolving...", total=len(queries))
        for q in queries:
            results.append(resolve(q, sources=src_list, query_type=query_type, unichem=not no_unichem, timeout=timeout))
            progress.advance(task)
    n_found = sum(1 for r in results if r.found)
    n_agree = sum(1 for r in results if r.agreement == "agree")
    n_dis = sum(1 for r in results if r.agreement == "disagree")
    console.print(f"\n[bold]Resolved {n_found}/{len(results)}[/bold]  ·  agree {n_agree}  ·  disagree {n_dis}  ·  sources {', '.join(src_list)}")
    for r in results[:5]:
        _print_one(r)
    if len(results) > 5:
        console.print(f"[dim]... {len(results) - 5} more; use --output / --json to save all.[/dim]")
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["query", "source", "status", "source_id", "name", "smiles", "inchikey", "formula", "molecular_weight", "url", "agreement", "consensus_inchikey", "error"])
            for r in results:
                for o in r.outcomes:
                    rec = r.by_source(o.source)
                    w.writerow([r.query, o.source, o.status, rec.source_id if rec else "", rec.name if rec else "", rec.smiles if rec else "",
                                rec.inchikey if rec else "", rec.formula if rec else "", rec.molecular_weight if rec else "", rec.url if rec else "",
                                r.agreement, r.consensus_inchikey or "", o.error or ""])
        console.print(f"[green]Saved:[/green] {output}")
    if json_out:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(json.dumps([r.model_dump() for r in results], indent=2, default=str))
        console.print(f"[green]Saved:[/green] {json_out}")


# =============================================================================
# init command
# =============================================================================

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

    console.print(table)


@app.command(name="lookup")
def lookup_cmd(
    query: str = typer.Argument(..., help="SMILES, PubChem CID, InChIKey, or Compound Name"),
    search_type: str = typer.Option("auto", "--type", "-t", help="Search type (auto, cid, smiles, name, inchikey)"),
    no_cache: bool = typer.Option(False, "--no-cache", help="Bypass local cache"),
) -> None:
    """Look up compound information by SMILES, CID, InChIKey, or Name."""
    _valid_types = {"auto", "smiles", "cid", "name", "inchikey"}
    if search_type.lower() not in _valid_types:
        console.print(f"[red]Error:[/red] Unknown --type {search_type!r}. Valid: {sorted(_valid_types)}")
        raise typer.Exit(code=1)
    with console.status(f"[bold green]Searching for '{query}'...[/bold green]"):
        compound = lookup(query, search_type=search_type, use_cache=not no_cache)

    if not compound:
        console.print(f"[red]✖ No compound found for query:[/red] [bold]{query}[/bold]")
        raise typer.Exit(code=1)

    table = Table(title=f"Compound Details: {query}", show_header=True, header_style="bold blue")
    table.add_column("Property", style="cyan", width=24)
    table.add_column("Value", style="white")

    if compound.cid:
        table.add_row("PubChem CID", str(compound.cid))
    if compound.iupac_name:
        table.add_row("IUPAC Name", compound.iupac_name)
    if compound.canonical_smiles:
        table.add_row("Canonical SMILES", compound.canonical_smiles)
    if compound.isomeric_smiles and compound.isomeric_smiles != compound.canonical_smiles:
        table.add_row("Isomeric SMILES", compound.isomeric_smiles)
    if compound.molecular_formula:
        table.add_row("Formula", compound.molecular_formula)
    if compound.molecular_weight:
        table.add_row("Molecular Weight", f"{compound.molecular_weight:.4f} g/mol")
    if compound.xlogp is not None:
        table.add_row("XLogP", f"{compound.xlogp:.2f}")
    if compound.hbond_donor_count is not None:
        table.add_row("H-Bond Donors", str(compound.hbond_donor_count))
    if compound.hbond_acceptor_count is not None:
        table.add_row("H-Bond Acceptors", str(compound.hbond_acceptor_count))
    if compound.inchi:
        table.add_row("InChI", compound.inchi)
    if compound.inchikey:
        table.add_row("InChIKey", compound.inchikey)

    console.print(table)


@app.command()
def batch(
    input_file: Path = typer.Argument(..., help="Path to input file (.csv, .tsv, .xlsx, .smi, .sdf)"),
    output_file: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output file path",
    ),
    format: str = typer.Option("csv", "--format", "-f", help="Output format (csv, xlsx, json)"),
    keep_duplicates: bool = typer.Option(False, "--keep-duplicates", help="Do not remove duplicate entries"),
) -> None:
    """Process a batch of SMILES from a file and retrieve metadata."""
    if not input_file.exists():
        console.print(f"[red]✖ Input file not found: {input_file}[/red]")
        raise typer.Exit(code=1)

    if not output_file:
        output_file = input_file.with_name(f"{input_file.stem}_results.{format}")

    report_file = output_file.with_name(f"{output_file.stem}_report.log")

    try:
        queries = parse_compounds_file(input_file)
        original_count = len(queries)

        if not keep_duplicates:
            queries = list(dict.fromkeys(queries))

        console.print(f"[blue]ℹ[/blue] Loaded {original_count} compounds ({len(queries)} unique).")
    except Exception as e:
        console.print(f"[red]✖ Failed to parse input file: {e}[/red]")
        raise typer.Exit(code=1)

    results = []
    log_entries = []
    found_count = 0

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        console=console,
    ) as progress:
        task = progress.add_task("[cyan]Looking up compounds...", total=len(queries))

        for query in queries:
            compound = lookup(query, use_cache=True)
            if compound:
                results.append(compound)
                found_count += 1
                log_entries.append(f"[SUCCESS] Query: '{query}' -> Found CID: {compound.cid}")
            else:
                results.append(PubChemCompound(input_query=query))
                log_entries.append(f"[FAILED]  Query: '{query}' -> Reason: Not found in PubChem")

            progress.advance(task)

    try:
        with open(report_file, "w", encoding="utf-8") as f:
            f.write("\n".join(log_entries))

        export_results(results, output_file, format)
        console.print("\n[green]✔[/green] Batch processing complete!")
        console.print(f"  • Found: [green]{found_count}[/green]")
        console.print(f"  • Missed: [red]{len(queries) - found_count}[/red]")
        console.print(f"  • Results saved to: [cyan]{output_file}[/cyan]")
        console.print(f"  • Report log saved to: [cyan]{report_file}[/cyan]")
    except Exception as e:
        console.print(f"[red]✖ Failed to export results: {e}[/red]")


def _sanitize_name_for_filename(name: str) -> str:
    """Generate a clean, filesystem-safe filename stem."""
    clean = re.sub(r'[\\/*?:"<>|]', '_', str(name).strip())
    clean = re.sub(r'_+', '_', clean).strip('_')
    if not clean:
        clean = "compound"
    if len(clean) > 40:
        short_hash = hashlib.md5(str(name).encode("utf-8")).hexdigest()[:8]
        clean = clean[:30] + "_" + short_hash
    return clean


def _resolve_smiles_for_query(query_str: str, offline: bool = False) -> tuple[Optional[str], Optional[int], Optional[str]]:
    """
    Resolve SMILES, CID, and compound title for a query string.
    Returns (smiles, cid, title).

    With ``offline=True`` a valid SMILES is accepted as-is and PubChem is never contacted;
    a name or CID still requires a lookup and returns (None, None, None) if unresolved offline.
    """
    val = validate_smiles(query_str)
    if val.is_valid and val.canonical_smiles:
        if offline:
            return val.canonical_smiles, None, val.canonical_smiles
        # Query itself is a valid SMILES string
        # Check if CID can also be found in cache/PubChem without failing if not found
        compound = lookup(val.canonical_smiles, use_cache=True)
        cid = compound.cid if compound else None
        title = (compound.iupac_name or val.canonical_smiles) if compound else val.canonical_smiles
        return val.canonical_smiles, cid, title
    if offline:
        return None, None, None

    # Not directly a SMILES; lookup in PubChem/cache
    compound = lookup(query_str, use_cache=True)
    if compound:
        smi = compound.canonical_smiles or compound.isomeric_smiles
        return smi, compound.cid, compound.iupac_name or query_str

    return None, None, None


@app.command()
def download(
    query: Optional[str] = typer.Argument(None, help="PubChem CID, SMILES, or compound name"),
    input_file: Optional[Path] = typer.Option(None, "--file", "-i", help="Batch download from file (CSVs, TSV, XLSX, SMI, SDF)"),
    format: str = typer.Option("sdf", "--format", "-f", help="Structure format (sdf, mol, pdb, png)"),
    is_3d: bool = typer.Option(True, "--3d/--2d", help="Download/generate 3D structure (use --2d for 2D)"),
    output_dir: Path = typer.Option(Path("structures"), "--output-dir", "-o", help="Output directory"),
    force: bool = typer.Option(False, "--force", help="Overwrite existing files (disables resume)"),
    gen: Optional[str] = typer.Option(
        None,
        "--gen",
        "-g",
        help="Structure generation mode: 'all' (generate all locally from SMILES) or 'missing' (generate only when not in PubChem/database)",
    ),
) -> None:
    """Download chemical structures from PubChem or generate 2D/3D structures locally from SMILES."""
    if not query and not input_file:
        console.print("[red]✖ You must provide either a query (CID/SMILES/Name) or an --file input.[/red]")
        raise typer.Exit(code=1)

    if gen is not None:
        gen = gen.lower().strip()
        if gen not in ["all", "missing"]:
            console.print(f"[red]✖ Invalid value for --gen: '{gen}'. Allowed values are 'all' or 'missing'.[/red]")
            raise typer.Exit(code=1)

    fmt = format.lower().strip()
    dimension = "3d" if is_3d else "2d"

    if gen == "all" and fmt not in ["sdf", "mol", "pdb"]:
        console.print(f"[red]✖ Structure generation does not support format '{format}'. Supported formats: sdf, mol, pdb.[/red]")
        raise typer.Exit(code=1)

    out_dir_str = str(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------------------
    # 1. Single Query Handling
    # -------------------------------------------------------------------------
    if query:
        query_str = str(query).strip()

        if gen == "all":
            with console.status(f"[bold green]Generating {dimension.upper()} {fmt.upper()} for '{query_str}'...[/bold green]"):
                smi, cid_val, title = _resolve_smiles_for_query(query_str, offline=True)
                if smi:
                    stem = str(cid_val) if cid_val else _sanitize_name_for_filename(query_str)
                    out_path = output_dir / f"{stem}_{dimension}.{fmt}"
                    status = generate_structure(smi, out_path, format=fmt, dimension=dimension, force=force, title=title)
                    if status == "Generated":
                        console.print(f"[green]✔[/green] Successfully generated to [cyan]{out_path}[/cyan]")
                    elif "Skipped" in status:
                        console.print(f"[blue]ℹ[/blue] {status}")
                    else:
                        console.print(f"[red]✖[/red] Failed to generate: {status}")
                        raise typer.Exit(code=1)
                else:
                    console.print(f"[red]✖[/red] Could not resolve SMILES for '{query_str}' to generate structure.")
                    raise typer.Exit(code=1)

        elif gen == "missing":
            with console.status(f"[bold green]Downloading/Generating {dimension.upper()} {fmt.upper()} for '{query_str}'...[/bold green]"):
                target_cid = int(query_str) if query_str.isdigit() else None
                compound = None
                if not target_cid:
                    compound = lookup(query_str, use_cache=True)
                    if compound and compound.cid:
                        target_cid = compound.cid

                downloaded = False
                if target_cid:
                    dl_status = download_structure(target_cid, fmt, dimension, out_dir_str, force)
                    if dl_status == "Downloaded":
                        console.print(f"[green]✔[/green] Successfully downloaded to [cyan]{output_dir}/{target_cid}_{dimension}.{fmt}[/cyan]")
                        downloaded = True
                    elif "Skipped" in dl_status:
                        console.print(f"[blue]ℹ[/blue] {dl_status}")
                        downloaded = True

                if not downloaded:
                    # Fallback to local generation from SMILES
                    smi = None
                    title = query_str
                    cid_val = target_cid

                    if compound and (compound.canonical_smiles or compound.isomeric_smiles):
                        smi = compound.canonical_smiles or compound.isomeric_smiles
                        title = compound.iupac_name or query_str
                    else:
                        val = validate_smiles(query_str)
                        if val.is_valid and val.canonical_smiles:
                            smi = val.canonical_smiles

                    if smi and fmt in ["sdf", "mol", "pdb"]:
                        stem = str(cid_val) if cid_val else _sanitize_name_for_filename(query_str)
                        out_path = output_dir / f"{stem}_{dimension}.{fmt}"
                        gen_status = generate_structure(smi, out_path, format=fmt, dimension=dimension, force=force, title=title)
                        if gen_status == "Generated":
                            console.print(f"[green]✔[/green] (PubChem unavailable) Generated structure locally to [cyan]{out_path}[/cyan]")
                        elif "Skipped" in gen_status:
                            console.print(f"[blue]ℹ[/blue] {gen_status}")
                        else:
                            console.print(f"[red]✖[/red] Failed to generate structure: {gen_status}")
                            raise typer.Exit(code=1)
                    else:
                        console.print("[red]✖[/red] Structure not found in PubChem/database and could not be generated locally.")
                        raise typer.Exit(code=1)

        else:
            # Default: PubChem download
            with console.status(f"[bold green]Downloading {dimension.upper()} {fmt.upper()} for '{query_str}'...[/bold green]"):
                target_cid = int(query_str) if query_str.isdigit() else None
                if not target_cid:
                    compound = lookup(query_str, use_cache=True)
                    if compound and compound.cid:
                        target_cid = compound.cid

                if target_cid:
                    status = download_structure(target_cid, fmt, dimension, out_dir_str, force)
                    if status == "Downloaded":
                        console.print(f"[green]✔[/green] Successfully downloaded to [cyan]{output_dir}/{target_cid}_{dimension}.{fmt}[/cyan]")
                    elif "Skipped" in status:
                        console.print(f"[blue]ℹ[/blue] {status}")
                    else:
                        console.print(f"[red]✖[/red] Failed to download: {status}")
                        raise typer.Exit(code=1)
                else:
                    console.print(f"[red]✖[/red] Could not resolve PubChem CID for '{query_str}'. Use --gen all or --gen missing to generate from SMILES.")
                    raise typer.Exit(code=1)

    # -------------------------------------------------------------------------
    # 2. Batch Download/Generation from File
    # -------------------------------------------------------------------------
    if input_file:
        if not input_file.exists():
            console.print(f"[red]✖ Input file not found: {input_file}[/red]")
            raise typer.Exit(code=1)

        report_file = output_dir / f"{input_file.stem}_download_report.log"

        try:
            queries = parse_compounds_file(input_file)
            queries = list(dict.fromkeys(queries))
        except Exception as e:
            console.print(f"[red]✖ Failed to parse input file: {e}[/red]")
            raise typer.Exit(code=1)

        action_desc = "generating" if gen == "all" else "downloading/generating" if gen == "missing" else "downloading"
        console.print(f"[blue]ℹ[/blue] Batch {action_desc} structures for {len(queries)} unique compounds...")

        downloaded_count = 0
        generated_count = 0
        skipped_count = 0
        failed_count = 0
        log_entries = []

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            console=console,
        ) as progress:
            task = progress.add_task("[cyan]Processing structures...", total=len(queries))

            for q in queries:
                q_str = str(q).strip()

                if gen == "all":
                    smi, cid_val, title = _resolve_smiles_for_query(q_str, offline=True)
                    if smi:
                        stem = str(cid_val) if cid_val else _sanitize_name_for_filename(q_str)
                        out_path = output_dir / f"{stem}_{dimension}.{fmt}"
                        status = generate_structure(smi, out_path, format=fmt, dimension=dimension, force=force, title=title)
                        if status == "Generated":
                            generated_count += 1
                            log_entries.append(f"[GENERATED] Query: '{q_str}' -> Saved to {out_path.name}")
                        elif "Skipped" in status:
                            skipped_count += 1
                            log_entries.append(f"[SKIPPED]   Query: '{q_str}' -> Reason: {status}")
                        else:
                            failed_count += 1
                            log_entries.append(f"[FAILED]    Query: '{q_str}' -> Reason: {status}")
                    else:
                        failed_count += 1
                        log_entries.append(f"[FAILED]    Query: '{q_str}' -> Reason: Could not resolve SMILES")

                elif gen == "missing":
                    target_cid = int(q_str) if q_str.isdigit() else None
                    compound = None
                    if not target_cid:
                        compound = lookup(q_str, use_cache=True)
                        if compound and compound.cid:
                            target_cid = compound.cid

                    handled = False
                    if target_cid:
                        dl_status = download_structure(target_cid, fmt, dimension, out_dir_str, force)
                        if dl_status == "Downloaded":
                            downloaded_count += 1
                            log_entries.append(f"[DOWNLOADED] Query: '{q_str}' -> CID {target_cid} -> Saved from PubChem")
                            handled = True
                        elif "Skipped" in dl_status:
                            skipped_count += 1
                            log_entries.append(f"[SKIPPED]    Query: '{q_str}' -> CID {target_cid} -> Reason: {dl_status}")
                            handled = True

                    if not handled:
                        # Fall back to generate locally from SMILES
                        smi = None
                        title = q_str
                        cid_val = target_cid
                        if compound and (compound.canonical_smiles or compound.isomeric_smiles):
                            smi = compound.canonical_smiles or compound.isomeric_smiles
                            title = compound.iupac_name or q_str
                        else:
                            val = validate_smiles(q_str)
                            if val.is_valid and val.canonical_smiles:
                                smi = val.canonical_smiles

                        if smi and fmt in ["sdf", "mol", "pdb"]:
                            stem = str(cid_val) if cid_val else _sanitize_name_for_filename(q_str)
                            out_path = output_dir / f"{stem}_{dimension}.{fmt}"
                            gen_status = generate_structure(smi, out_path, format=fmt, dimension=dimension, force=force, title=title)
                            if gen_status == "Generated":
                                generated_count += 1
                                log_entries.append(f"[GENERATED]  Query: '{q_str}' -> (PubChem unavailable) Generated locally to {out_path.name}")
                            elif "Skipped" in gen_status:
                                skipped_count += 1
                                log_entries.append(f"[SKIPPED]    Query: '{q_str}' -> Reason: {gen_status}")
                            else:
                                failed_count += 1
                                log_entries.append(f"[FAILED]     Query: '{q_str}' -> Reason: {gen_status}")
                        else:
                            failed_count += 1
                            log_entries.append(f"[FAILED]     Query: '{q_str}' -> Reason: Not available in PubChem and cannot generate")

                else:
                    # Default: PubChem download only
                    target_cid = int(q_str) if q_str.isdigit() else None
                    if not target_cid:
                        compound = lookup(q_str, use_cache=True)
                        if compound and compound.cid:
                            target_cid = compound.cid

                    if target_cid:
                        status = download_structure(target_cid, fmt, dimension, out_dir_str, force)
                        if status == "Downloaded":
                            downloaded_count += 1
                            log_entries.append(f"[DOWNLOADED] Query: '{q_str}' -> CID {target_cid} -> Downloaded")
                        elif "Skipped" in status:
                            skipped_count += 1
                            log_entries.append(f"[SKIPPED]    Query: '{q_str}' -> Reason: {status}")
                        else:
                            failed_count += 1
                            log_entries.append(f"[FAILED]     Query: '{q_str}' -> Reason: {status}")
                    else:
                        failed_count += 1
                        log_entries.append(f"[FAILED]     Query: '{q_str}' -> Reason: Could not resolve PubChem CID")

                progress.advance(task)

        # Save status log
        with open(report_file, "w", encoding="utf-8") as f:
            f.write("\n".join(log_entries))

        console.print("\n[green]✔[/green] Batch processing complete!")
        if gen in ["all", "missing"]:
            console.print(f"    Downloaded: [green]{downloaded_count}[/green]")
            console.print(f"    Generated:  [green]{generated_count}[/green]")
        else:
            console.print(f"    Downloaded: [green]{downloaded_count}[/green]")
        console.print(f"    Skipped (Already exist): [blue]{skipped_count}[/blue]")
        console.print(f"    Failed/Not Found:        [red]{failed_count}[/red]")
        console.print(f"    Saved in: [cyan]{output_dir}[/cyan]")
        console.print(f"    Report log saved to: [cyan]{report_file}[/cyan]")



# ═══════════════════════════════════════════════════════════════════════════════
# fingerprint command
# ═══════════════════════════════════════════════════════════════════════════════

@app.command(name="fingerprint")
def fingerprint_cmd(
    smiles: str = typer.Argument(None, help="SMILES string of the compound to fingerprint"),
    file: Optional[Path] = typer.Option(None, "--file", "-i", help="Input file (.csv, .smi, .txt) for batch processing"),
    fp_type: str = typer.Option("ecfp4", "--type", "-t", help="Fingerprint type: ecfp4, ecfp6, fcfp4, maccs, rdkit, atompair, torsion, all"),
    bits: int = typer.Option(2048, "--bits", "-b", help="Number of bits (ignored for MACCS which is fixed at 167)"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output CSV file for batch results"),
) -> None:
    """Generate molecular fingerprints from SMILES (offline, RDKit-based).

    Supports ECFP4, ECFP6, FCFP4, MACCS, RDKit, AtomPair, and Topological Torsion fingerprints.
    """
    if smiles is None and file is None:
        console.print("[red]Error:[/red] Provide a SMILES argument or --file for batch input.")
        raise typer.Exit(code=1)

    valid_types = {"ecfp4", "ecfp6", "fcfp4", "maccs", "rdkit", "atompair", "torsion", "all"}
    if fp_type.lower() not in valid_types:
        console.print(f"[red]Error:[/red] Unknown fingerprint type '{fp_type}'. Valid: {sorted(valid_types)}")
        raise typer.Exit(code=1)

    # ── Single SMILES mode ──────────────────────────────────────────────────
    if smiles is not None and file is None:
        try:
            result = compute_fingerprint(smiles, fp_type=fp_type, n_bits=bits)
        except ValueError as e:
            console.print(f"[red]Error:[/red] {e}")
            raise typer.Exit(code=1)

        results = result if isinstance(result, list) else [result]

        table = Table(
            title=f"[bold]Fingerprint Results[/bold] — {smiles[:60]}{'...' if len(smiles) > 60 else ''}",
            show_header=True,
            header_style="bold magenta",
        )
        table.add_column("Type", style="cyan", min_width=12)
        table.add_column("Bits", justify="right", style="white")
        table.add_column("On Bits", justify="right", style="green")
        table.add_column("Density", justify="right", style="yellow")
        table.add_column("Bit String (preview)", style="dim")

        for r in results:
            preview = r.bit_string[:64] + "..." if len(r.bit_string) > 64 else r.bit_string
            table.add_row(
                r.fingerprint_type.upper(),
                str(r.n_bits),
                str(r.n_on_bits),
                f"{r.density:.4f}",
                preview,
            )
        console.print(table)
        if len(results) == 1:
            console.print(f"\n[dim]Full bit string ({results[0].n_bits} bits):[/dim]")
            console.print(f"[white]{results[0].bit_string}[/white]")
        return

    # ── Batch file mode ─────────────────────────────────────────────────────
    from chemlitmus.utils.parsers import parse_compounds_file
    queries = parse_compounds_file(file)

    rows = []
    errors = 0
    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
                  BarColumn(), TaskProgressColumn(), console=console) as progress:
        task = progress.add_task("[green]Computing fingerprints...", total=len(queries))
        for q in queries:
            try:
                res = compute_fingerprint(q.strip(), fp_type=fp_type, n_bits=bits)
                res_list = res if isinstance(res, list) else [res]
                for r in res_list:
                    rows.append({
                        "smiles": r.smiles,
                        "fp_type": r.fingerprint_type,
                        "n_bits": r.n_bits,
                        "n_on_bits": r.n_on_bits,
                        "density": r.density,
                        "bit_string": r.bit_string,
                        "hex_string": r.hex_string,
                    })
            except Exception:
                errors += 1
            progress.advance(task)

    if output:
        import csv
        output.parent.mkdir(parents=True, exist_ok=True)
        if rows:
            with open(output, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
        console.print(f"\n[green]Fingerprints saved to[/green]: [cyan]{output}[/cyan]")
    else:
        table = Table(show_header=True, header_style="bold magenta", title="Batch Fingerprint Results")
        table.add_column("SMILES", style="cyan", max_width=30, overflow="fold")
        table.add_column("Type", style="white")
        table.add_column("Bits", justify="right")
        table.add_column("On Bits", justify="right", style="green")
        table.add_column("Density", justify="right", style="yellow")
        for row in rows[:50]:
            table.add_row(row["smiles"][:28], row["fp_type"].upper(), str(row["n_bits"]), str(row["n_on_bits"]), f"{row['density']:.4f}")
        console.print(table)

    console.print(f"\n[green]Done.[/green] Processed: [white]{len(queries)}[/white] | Errors: [red]{errors}[/red]")


# ═══════════════════════════════════════════════════════════════════════════════
# similar command
# ═══════════════════════════════════════════════════════════════════════════════

@app.command(name="similar")
def similar_cmd(
    query: str = typer.Argument(..., help="Query SMILES string to search for"),
    file: Optional[Path] = typer.Option(None, "--file", "-i", help="Library file (.csv, .smi, .txt) to search against"),
    threshold: float = typer.Option(0.5, "--threshold", "-t", help="Minimum Tanimoto similarity (0.0-1.0)"),
    top: int = typer.Option(10, "--top", "-n", help="Number of top results to return"),
    fp_type: str = typer.Option("ecfp4", "--fp-type", help="Fingerprint type: ecfp4, ecfp6, fcfp4, maccs, rdkit, atompair, torsion"),
    bits: int = typer.Option(2048, "--bits", "-b", help="Number of fingerprint bits"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Save results to CSV file"),
) -> None:
    """Search a library file for compounds similar to a query SMILES (Tanimoto similarity).

    Ranks library compounds by Tanimoto similarity. Uses RDKit fingerprints offline.
    """
    if file is None:
        console.print("[red]Error:[/red] --file / -i is required for similarity search.")
        raise typer.Exit(code=1)

    valid_types = {"ecfp4", "ecfp6", "fcfp4", "maccs", "rdkit", "atompair", "torsion"}
    if fp_type.lower() not in valid_types:
        console.print(f"[red]Error:[/red] Unknown fp-type '{fp_type}'. Valid: {sorted(valid_types)}")
        raise typer.Exit(code=1)

    from chemlitmus.utils.parsers import parse_compounds_file
    library = parse_compounds_file(file)

    with console.status(f"[bold green]Searching {len(library)} library compounds...[/bold green]"):
        try:
            hits = compute_similarity(
                query,
                library,
                fp_type=fp_type.lower(),
                n_bits=bits,
                threshold=threshold,
                top_n=top,
            )
        except ValueError as e:
            console.print(f"[red]Error:[/red] {e}")
            raise typer.Exit(code=1)

    if not hits:
        console.print(f"[yellow]No compounds found with Tanimoto similarity >= {threshold}[/yellow]")
        raise typer.Exit()

    table = Table(
        title=f"[bold]Similarity Search Results[/bold] — Top {len(hits)} hits (threshold={threshold}, fp={fp_type.upper()})",
        show_header=True,
        header_style="bold magenta",
    )
    table.add_column("Rank", justify="right", style="dim", width=5)
    table.add_column("Similarity", justify="right", style="bold green", width=12)
    table.add_column("Hit SMILES", style="cyan")

    for hit in hits:
        sim_bar = "#" * int(hit.similarity * 10) + "-" * (10 - int(hit.similarity * 10))
        table.add_row(str(hit.rank), f"{hit.similarity:.4f} {sim_bar}", hit.hit)

    console.print(table)
    console.print(f"\n[dim]Query:[/dim] [white]{query}[/white]")
    console.print(f"[dim]Library:[/dim] [white]{file}[/white] ({len(library)} compounds)")
    console.print(f"[dim]Fingerprint:[/dim] [white]{fp_type.upper()}, {bits} bits[/white]")

    if output:
        import csv
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["rank", "query", "hit", "similarity", "fingerprint_type"])
            writer.writeheader()
            for h in hits:
                writer.writerow(h.model_dump())
        console.print(f"\n[green]Results saved to[/green]: [cyan]{output}[/cyan]")


# ═══════════════════════════════════════════════════════════════════════════════
# filter command
# ═══════════════════════════════════════════════════════════════════════════════

@app.command(name="filter")
def filter_cmd(
    smiles: str = typer.Argument(None, help="Single SMILES string to evaluate"),
    file: Optional[Path] = typer.Option(None, "--file", "-i", help="Input file for batch evaluation"),
    rules: str = typer.Option("all", "--rules", "-r", help="Comma-separated rules: lipinski,veber,ghose,egan,ro3,pains,qed or all"),
    fail: bool = typer.Option(False, "--fail", help="Invert output: keep only compounds that FAIL the filter (useful for PAINS removal)"),
    qed_min: float = typer.Option(0.0, "--qed-min", help="Minimum QED score to keep (0.0-1.0)"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Save results to CSV file"),
    prep: str = typer.Option("implicit-h", "--prep", help="Molecule preparation for substructure-based rules: implicit-h (default), explicit-h or kekule. Recorded in the output."),
) -> None:
    """Apply drug-likeness and ADMET filters to SMILES compounds (offline, RDKit-based).

    Evaluates: Lipinski Ro5, Veber, Ghose, Egan, Ro3, PAINS alerts, and QED score.
    Substructure alerts depend on how molecules are prepared; --prep makes that choice
    explicit and writes it into the results so a screen can be reproduced.
    """
    if smiles is None and file is None:
        console.print("[red]Error:[/red] Provide a SMILES argument or --file for batch input.")
        raise typer.Exit(code=1)

    # Parse rules list
    rule_list = [r.strip().lower() for r in rules.split(",")]

    # ── Single SMILES mode ──────────────────────────────────────────────────
    if smiles is not None and file is None:
        try:
            result = apply_filters(smiles, rules=rule_list, preparation=prep)
        except ValueError as e:
            console.print(f"[red]Error:[/red] {e}")
            raise typer.Exit(code=1)

        if result.error:
            console.print(f"[red]Error:[/red] {result.error}")
            raise typer.Exit(code=1)

        console.print(f"[dim]Preparation for substructure rules: {result.preparation}[/dim]")

        # Properties table
        prop_table = Table(title="[bold]Physicochemical Properties[/bold]", show_header=True, header_style="bold blue")
        prop_table.add_column("Property", style="cyan")
        prop_table.add_column("Value", justify="right", style="white")
        prop_table.add_column("Property", style="cyan")
        prop_table.add_column("Value", justify="right", style="white")

        prop_table.add_row(
            "Mol Weight (g/mol)", f"{result.molecular_weight:.3f}",
            "LogP (MolLogP)", f"{result.logp:.3f}",
        )
        prop_table.add_row(
            "H-Bond Donors", str(result.hbd),
            "H-Bond Acceptors", str(result.hba),
        )
        prop_table.add_row(
            "TPSA (A^2)", f"{result.tpsa:.3f}",
            "Rotatable Bonds", str(result.rotatable_bonds),
        )
        prop_table.add_row(
            "Heavy Atoms", str(result.heavy_atom_count),
            "Molar Refractivity", f"{result.molar_refractivity:.3f}",
        )
        prop_table.add_row(
            "QED Score", f"[{'green' if result.qed_score >= 0.5 else 'yellow'}]{result.qed_score:.4f}[/{'green' if result.qed_score >= 0.5 else 'yellow'}]",
            "", "",
        )
        console.print(prop_table)

        # Filter results table
        filt_table = Table(title="[bold]Filter Results[/bold]", show_header=True, header_style="bold magenta")
        filt_table.add_column("Rule", style="cyan", min_width=18)
        filt_table.add_column("Result", width=10)
        filt_table.add_column("Details", style="dim")

        rule_map = {
            "lipinski": ("Lipinski Ro5", result.lipinski),
            "veber": ("Veber", result.veber),
            "ghose": ("Ghose", result.ghose),
            "egan": ("Egan", result.egan),
            "ro3": ("Ro3 (Lead-like)", result.ro3),
            "pains": ("PAINS", result.pains),
        }
        for key, (label, rule_result) in rule_map.items():
            if rule_result is not None:
                icon = "[green]PASS[/green]" if rule_result.passed else "[red]FAIL[/red]"
                filt_table.add_row(label, icon, rule_result.details)

        overall = "[green]PASS[/green]" if result.passes_all else "[red]FAIL[/red]"
        filt_table.add_row("[bold]Overall[/bold]", overall, "All requested rules")
        console.print(filt_table)
        return

    # ── Batch file mode ─────────────────────────────────────────────────────
    import csv
    from chemlitmus.utils.parsers import parse_compounds_file
    queries = parse_compounds_file(file)

    rows = []
    pass_count = 0
    fail_count = 0
    error_count = 0

    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"),
                  BarColumn(), TaskProgressColumn(), console=console) as progress:
        task = progress.add_task("[green]Applying filters...", total=len(queries))
        for q in queries:
            try:
                r = apply_filters(q.strip(), rules=rule_list, preparation=prep)
                if r.error:
                    error_count += 1
                    progress.advance(task)
                    continue

                # Apply QED threshold
                if qed_min > 0 and (r.qed_score is None or r.qed_score < qed_min):
                    progress.advance(task)
                    continue

                is_kept = (r.passes_all and not fail) or (not r.passes_all and fail)
                if is_kept:
                    row = {
                        "smiles": r.smiles,
                        "mw": r.molecular_weight,
                        "logp": r.logp,
                        "hbd": r.hbd,
                        "hba": r.hba,
                        "tpsa": r.tpsa,
                        "rotatable_bonds": r.rotatable_bonds,
                        "heavy_atoms": r.heavy_atom_count,
                        "qed": r.qed_score,
                        "lipinski": r.lipinski.passed if r.lipinski else "",
                        "veber": r.veber.passed if r.veber else "",
                        "ghose": r.ghose.passed if r.ghose else "",
                        "egan": r.egan.passed if r.egan else "",
                        "ro3": r.ro3.passed if r.ro3 else "",
                        "pains_clean": r.pains.passed if r.pains else "",
                        "passes_all": r.passes_all,
                        "preparation": r.preparation,
                    }
                    rows.append(row)

                if r.passes_all:
                    pass_count += 1
                else:
                    fail_count += 1
            except ValueError as exc:
                # a bad --prep / --rules value is a usage error, not a per-compound failure
                progress.stop()
                console.print(f"[red]Error:[/red] {exc}")
                raise typer.Exit(code=1)
            except Exception as exc:
                error_count += 1
                logger.warning("filter failed for %r: %s", q, exc)
            progress.advance(task)

    # Show summary table (first 20)
    summary_table = Table(show_header=True, header_style="bold magenta",
                          title=f"[bold]Filter Results — {'Failures' if fail else 'Passes'} ({len(rows)} compounds)[/bold]")
    summary_table.add_column("SMILES", style="cyan", max_width=30, overflow="fold")
    summary_table.add_column("MW", justify="right")
    summary_table.add_column("LogP", justify="right")
    summary_table.add_column("QED", justify="right")
    summary_table.add_column("Overall", justify="center")

    for row in rows[:20]:
        overall_str = "[green]PASS[/green]" if row["passes_all"] else "[red]FAIL[/red]"
        summary_table.add_row(
            row["smiles"][:28],
            str(row["mw"]),
            str(row["logp"]),
            str(row["qed"]),
            overall_str,
        )
    console.print(summary_table)

    console.print(f"\n[bold]Summary[/bold]: Total={len(queries)} | Pass={pass_count} | Fail={fail_count} | Errors={error_count}")

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        if rows:
            with open(output, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
        console.print(f"[green]Results saved to[/green]: [cyan]{output}[/cyan]")





# =============================================================================
# standardize command
# =============================================================================

@app.command(name="standardize")
def standardize_cmd(
    smiles: str = typer.Argument(None, help="Input SMILES string to standardize."),
    steps: str = typer.Option(
        "all",
        "--steps", "-s",
        help=(
            "Comma-separated list of steps to apply. "
            "Valid: fragment, neutralize, tautomer, canonical, all. "
            "Default: all"
        ),
    ),
    file: Optional[Path] = typer.Option(None, "--file", "-f", help="Input file (CSV or .smi, one SMILES per line / column)."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Save results to CSV."),
    show_diff: bool = typer.Option(False, "--show-diff", help="Show per-step diff of SMILES changes."),
) -> None:
    """Standardize SMILES strings: salt stripping, neutralization, tautomer canonicalization.

    Uses RDKit MolStandardize — fully offline, no network access required.

    Steps applied (in order):

    \b
      fragment   - Keep largest organic fragment (salt stripping)
      neutralize - Neutralize charged atoms (e.g. carboxylate -> acid)
      tautomer   - Canonicalize tautomers to a single stable form
      canonical  - Generate canonical SMILES representation
    """
    import csv

    if smiles is None and file is None:
        console.print("[red]Error:[/red] Provide either a SMILES argument or --file.")
        raise typer.Exit(code=1)

    # Resolve steps list
    raw_steps = [s.strip().lower() for s in steps.split(",")]
    try:
        if raw_steps == ["all"]:
            step_list = None  # means all
        else:
            invalid = [s for s in raw_steps if s not in STANDARDIZE_STEPS]
            if invalid:
                console.print(f"[red]Invalid steps:[/red] {invalid}. Valid: {STANDARDIZE_STEPS + ['all']}")
                raise typer.Exit(code=1)
            step_list = raw_steps
    except typer.Exit:
        raise
    except Exception as exc:
        console.print(f"[red]Error:[/red] {exc}")
        raise typer.Exit(code=1)

    # Single SMILES mode
    if smiles and file is None:
        result = standardize_smiles(smiles, steps=step_list)
        _print_standardize_result(result, show_diff=show_diff)
        return

    # Batch file mode
    from chemlitmus.utils.parsers import parse_compounds_file
    compounds = parse_compounds_file(file)
    smiles_list = [c if isinstance(c, str) else c.get("smiles", "") for c in compounds]

    results = []
    with console.status(f"[bold green]Standardizing {len(smiles_list)} compounds...[/bold green]"):
        for smi in smiles_list:
            if not smi.strip():
                continue
            r = standardize_smiles(smi, steps=step_list)
            results.append(r)

    # Print summary table
    table = Table(title=f"Standardization Results ({len(results)} compounds)", show_lines=False)
    table.add_column("Input SMILES", style="dim", max_width=35)
    table.add_column("Output SMILES", style="bright_cyan", max_width=35)
    table.add_column("Changed", justify="center")
    table.add_column("Error", style="red", max_width=25)

    for r in results:
        changed_str = "[bright_green]Yes[/bright_green]" if r.changed else "[dim]No[/dim]"
        table.add_row(
            r.input_smiles[:35],
            (r.output_smiles or "")[:35],
            changed_str,
            r.error or "",
        )
    console.print(table)

    n_changed = sum(1 for r in results if r.changed)
    n_errors = sum(1 for r in results if r.error)
    console.print(
        f"\n[bold]Summary:[/bold] {len(results)} processed, "
        f"[bright_green]{n_changed} changed[/bright_green], "
        f"[red]{n_errors} errors[/red]"
    )

    if output:
        with open(output, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["input_smiles", "output_smiles", "changed", "steps_applied", "error"],
            )
            writer.writeheader()
            for r in results:
                writer.writerow({
                    "input_smiles": r.input_smiles,
                    "output_smiles": r.output_smiles or "",
                    "changed": r.changed,
                    "steps_applied": ",".join(r.steps_applied),
                    "error": r.error or "",
                })
        console.print(f"[green]Results saved to:[/green] {output}")


def _print_standardize_result(result: StandardizeResult, show_diff: bool = False) -> None:
    """Pretty-print a single StandardizeResult."""
    from rich.text import Text

    if result.error and not result.output_smiles:
        console.print(f"[red]Error:[/red] {result.error}")
        return

    # Header panel
    from rich.panel import Panel
    body = Text()
    body.append("\n")
    body.append("  Input:  ", style="dim white")
    body.append(result.input_smiles, style="white")
    body.append("\n")
    body.append("  Output: ", style="dim white")
    body.append(result.output_smiles or "N/A", style="bold bright_cyan")
    body.append("\n")
    changed_style = "bold bright_green" if result.changed else "dim"
    changed_label = "Yes - structure was modified" if result.changed else "No - already canonical"
    body.append("  Changed: ", style="dim white")
    body.append(changed_label, style=changed_style)
    body.append("\n")
    if result.error:
        body.append("  Warning: ", style="yellow")
        body.append(result.error, style="yellow")
        body.append("\n")

    border = "bright_green" if not result.error else "yellow"
    console.print(Panel(body, title="[bold bright_cyan]Standardization Result[/bold bright_cyan]",
                         border_style=border, padding=(0, 2)))

    if show_diff and result.step_results:
        console.print()
        table = Table(title="Step-by-Step Changes", show_lines=True)
        table.add_column("Step", style="bold white", width=12)
        table.add_column("Input", style="dim", max_width=40)
        table.add_column("Output", style="bright_cyan", max_width=40)
        table.add_column("Changed", justify="center", width=8)
        for sr in result.step_results:
            changed_cell = "[bright_green]Yes[/bright_green]" if sr.changed else "[dim]No[/dim]"
            table.add_row(sr.step.capitalize(), sr.input_smiles, sr.output_smiles, changed_cell)
        console.print(table)


# =============================================================================
# iupacname command
# =============================================================================

@app.command(name="iupacname")
def iupacname_cmd(
    smiles: str = typer.Argument(None, help="Input SMILES string."),
    online: bool = typer.Option(
        False, "--online",
        help="Fetch IUPAC preferred name from PubChem API (cached for future offline use).",
    ),
    file: Optional[Path] = typer.Option(None, "--file", "-f", help="Input file (CSV or .smi, one SMILES per line)."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Save results to CSV."),
) -> None:
    """Generate IUPAC identifiers from SMILES strings (offline: InChI/InChIKey; optional: IUPAC name via PubChem).

    \b
    Offline output (always available):
      - Canonical SMILES
      - Molecular Formula
      - Molecular Weight
      - InChI string
      - InChIKey

    IUPAC name (requires --online on first use, then cached offline):
      Fetches the IUPAC preferred name from PubChem and stores it in the
      local SQLite cache. Subsequent calls for the same compound are offline.
    """
    import csv

    if smiles is None and file is None:
        console.print("[red]Error:[/red] Provide either a SMILES argument or --file.")
        raise typer.Exit(code=1)

    # Single SMILES mode
    if smiles and file is None:
        with console.status("[bold green]Computing identifiers...[/bold green]"):
            result = get_iupac_name(smiles, use_online=online)
        _print_iupacname_result(result)
        return

    # Batch mode
    from chemlitmus.utils.parsers import parse_compounds_file
    compounds = parse_compounds_file(file)
    smiles_list = [c if isinstance(c, str) else c.get("smiles", "") for c in compounds]

    results = []
    status_label = "Fetching identifiers (+ PubChem names)..." if online else "Computing InChI identifiers..."
    with console.status(f"[bold green]{status_label}[/bold green]"):
        for smi in smiles_list:
            if not smi.strip():
                continue
            r = get_iupac_name(smi, use_online=online)
            results.append(r)

    # Summary table
    table = Table(title=f"IUPAC Identifier Results ({len(results)} compounds)", show_lines=False)
    table.add_column("Input SMILES", style="dim", max_width=28)
    table.add_column("Formula", style="white", width=12)
    table.add_column("InChIKey", style="dim", width=28)
    table.add_column("IUPAC Name", style="bright_cyan", max_width=35)
    table.add_column("Source", style="dim", width=12)

    for r in results:
        if r.error:
            table.add_row(r.input_smiles[:28], "[red]ERROR[/red]", "", r.error[:35], "")
        else:
            table.add_row(
                r.input_smiles[:28],
                r.molecular_formula or "",
                r.inchikey or "",
                r.iupac_name or "[dim]—[/dim]",
                r.iupac_name_source,
            )
    console.print(table)

    if online:
        n_named = sum(1 for r in results if r.iupac_name)
        console.print(f"\n[dim]IUPAC names found: {n_named}/{len(results)}[/dim]")

    if output:
        with open(output, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["input_smiles", "canonical_smiles", "molecular_formula",
                            "molecular_weight", "inchi", "inchikey", "iupac_name", "iupac_name_source", "error"],
            )
            writer.writeheader()
            for r in results:
                writer.writerow({
                    "input_smiles": r.input_smiles,
                    "canonical_smiles": r.canonical_smiles or "",
                    "molecular_formula": r.molecular_formula or "",
                    "molecular_weight": r.molecular_weight or "",
                    "inchi": r.inchi or "",
                    "inchikey": r.inchikey or "",
                    "iupac_name": r.iupac_name or "",
                    "iupac_name_source": r.iupac_name_source,
                    "error": r.error or "",
                })
        console.print(f"[green]Results saved to:[/green] {output}")


def _print_iupacname_result(result: IUPACResult) -> None:
    """Pretty-print a single IUPACResult."""
    from rich.panel import Panel
    from rich.text import Text

    if result.error and not result.canonical_smiles:
        console.print(f"[red]Error:[/red] {result.error}")
        return

    body = Text()
    body.append("\n")
    body.append("  Input SMILES:    ", style="white")
    body.append(result.input_smiles, style="dim")
    body.append("\n")
    body.append("  Canonical SMILES:", style="white")
    body.append(f" {result.canonical_smiles}", style="bright_cyan")
    body.append("\n")
    body.append("  Formula:         ", style="white")
    body.append(result.molecular_formula or "N/A", style="bold white")
    body.append("\n")
    body.append("  MW (exact):      ", style="white")
    body.append(f"{result.molecular_weight} g/mol" if result.molecular_weight else "N/A", style="white")
    body.append("\n")
    body.append("\n")
    body.append("  InChI:           ", style="dim white")
    body.append((result.inchi or "N/A")[:65], style="dim")
    if result.inchi and len(result.inchi) > 65:
        body.append("...", style="dim")
    body.append("\n")
    body.append("  InChIKey:        ", style="dim white")
    body.append(result.inchikey or "N/A", style="bright_blue")
    body.append("\n")
    body.append("\n")
    body.append("  IUPAC Name:      ", style="white")
    if result.iupac_name:
        body.append(result.iupac_name, style="bold bright_green")
        body.append(f"  (source: {result.iupac_name_source})", style="dim")
    else:
        body.append("Not available offline. Use --online to fetch from PubChem.", style="dim yellow")
    body.append("\n")

    console.print(Panel(body, title="[bold bright_cyan]IUPAC Identifier Result[/bold bright_cyan]",
                         border_style="bright_cyan", padding=(0, 2)))




# =============================================================================
# substructure command
# =============================================================================

@app.command(name="substructure")
def substructure_cmd(
    query: str = typer.Argument(..., help="The substructure query (SMARTS or SMILES)."),
    file: Path = typer.Option(..., "--file", "-f", help="Library file to search (CSV or .smi, one SMILES per line)."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Save matching SMILES to CSV."),
    smiles_query: bool = typer.Option(False, "--smiles-query", help="Treat the query strictly as a SMILES string (default assumes SMARTS)."),
    prep: str = typer.Option("implicit-h", "--prep", help="Library molecule preparation: implicit-h (default), explicit-h or kekule. Recorded in the output."),
) -> None:
    """Search a library for compounds containing a specific substructure.
    
    Uses RDKit's HasSubstructMatch. Query is assumed to be SMARTS by default 
    for maximum expressive power (e.g. `[#9,#17]` for halogen).
    """
    import csv
    from chemlitmus.utils.parsers import parse_compounds_file

    compounds = parse_compounds_file(file)
    library_smiles = [c if isinstance(c, str) else c.get("smiles", "") for c in compounds]
    
    with console.status(f"[bold green]Searching {len(library_smiles)} compounds for substructure...[/bold green]"):
        try:
            hits = substructure_search(query=query, library=library_smiles, is_smarts=not smiles_query, preparation=prep)
        except Exception as exc:
            console.print(f"[red]Error:[/red] {exc}")
            raise typer.Exit(code=1)

    table = Table(title=f"Substructure Search Results  [dim](preparation: {prep})[/dim]", show_lines=False)
    table.add_column("Match SMILES", style="bright_cyan")
    table.add_column("Matched Atoms (indices)", style="dim")
    
    for hit in hits:
        table.add_row(
            hit.smiles,
            str(hit.match_indices)
        )
    console.print(table)
    
    query_type = "SMILES" if smiles_query else "SMARTS"
    console.print(
        f"\n[bold]Summary:[/bold] Found [bright_green]{len(hits)}[/bright_green] matches out of {len(library_smiles)} "
        f"compounds for {query_type} query: '{query}'."
    )
    
    if output:
        with open(output, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["smiles", "match_indices", "preparation"],
            )
            writer.writeheader()
            for hit in hits:
                writer.writerow({
                    "smiles": hit.smiles,
                    "match_indices": ",".join(map(str, hit.match_indices)),
                    "preparation": hit.preparation,
                })
        console.print(f"[green]Results saved to:[/green] {output}")

# =============================================================================
# tautomers command
# =============================================================================

@app.command(name="tautomers")
def tautomers_cmd(
    smiles: str = typer.Argument(None, help="Input SMILES string."),
    max_tautomers: int = typer.Option(1000, "--max", "-m", help="Maximum number of tautomers to generate."),
    file: Optional[Path] = typer.Option(None, "--file", "-f", help="Input file (CSV or .smi, one SMILES per line)."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Save results to CSV (one row per tautomer)."),
) -> None:
    """Enumerate all plausible tautomers for a given SMILES string.
    
    Critical for protein-ligand docking preparation as different tautomers 
    can bind with vastly different affinities. Fully offline via RDKit.
    """
    import csv

    if smiles is None and file is None:
        console.print("[red]Error:[/red] Provide either a SMILES argument or --file.")
        raise typer.Exit(code=1)

    # Single SMILES mode
    if smiles and file is None:
        result = enumerate_tautomers(smiles, max_tautomers=max_tautomers)
        _print_tautomer_result(result)
        if result.error:
            raise typer.Exit(code=1)
        return

    # Batch file mode
    from chemlitmus.utils.parsers import parse_compounds_file
    compounds = parse_compounds_file(file)
    smiles_list = [c if isinstance(c, str) else c.get("smiles", "") for c in compounds]

    results = []
    with console.status(f"[bold green]Enumerating tautomers for {len(smiles_list)} compounds...[/bold green]"):
        for smi in smiles_list:
            if not smi.strip():
                continue
            r = enumerate_tautomers(smi, max_tautomers=max_tautomers)
            results.append(r)

    # Print summary table
    table = Table(title=f"Tautomer Enumeration Results ({len(results)} compounds)", show_lines=False)
    table.add_column("Input SMILES", style="dim", max_width=35)
    table.add_column("Tautomers Found", justify="right", style="bright_cyan")
    table.add_column("Error", style="red", max_width=25)

    total_tauts = 0
    for r in results:
        table.add_row(
            r.input_smiles[:35],
            str(r.num_tautomers) if r.success else "-",
            r.error or "",
        )
        total_tauts += r.num_tautomers
    console.print(table)

    console.print(
        f"\n[bold]Summary:[/bold] {len(results)} inputs generated "
        f"[bright_green]{total_tauts} total tautomers[/bright_green]."
    )

    if output:
        # Explode mode: one row per tautomer
        with open(output, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["input_smiles", "tautomer_smiles", "is_canonical", "error"],
            )
            writer.writeheader()
            for r in results:
                if not r.success:
                    writer.writerow({
                        "input_smiles": r.input_smiles,
                        "tautomer_smiles": "",
                        "is_canonical": "",
                        "error": r.error,
                    })
                else:
                    for t_smi in r.tautomers:
                        writer.writerow({
                            "input_smiles": r.input_smiles,
                            "tautomer_smiles": t_smi,
                            "is_canonical": str(t_smi == r.canonical_tautomer),
                            "error": "",
                        })
        console.print(f"[green]Results saved to:[/green] {output} [dim](exploded to {total_tauts} rows)[/dim]")


def _print_tautomer_result(result: TautomerResult) -> None:
    """Pretty-print a single TautomerResult."""
    from rich.panel import Panel
    from rich.text import Text

    if not result.success:
        console.print(f"[red]Error:[/red] {result.error}")
        return

    body = Text()
    body.append("\n")
    body.append("  Input:           ", style="dim white")
    body.append(result.input_smiles, style="white")
    body.append("\n")
    body.append("  Tautomers found: ", style="dim white")
    body.append(str(result.num_tautomers), style="bold bright_cyan")
    body.append("\n")

    console.print(Panel(body, title="[bold bright_cyan]Tautomer Enumeration[/bold bright_cyan]",
                         border_style="bright_cyan", padding=(0, 2)))

    if result.tautomers:
        console.print()
        table = Table(title=f"All {result.num_tautomers} Tautomers", show_lines=True)
        table.add_column("#", style="dim", justify="right")
        table.add_column("SMILES", style="white")
        table.add_column("Canonical", justify="center")
        
        for i, t_smi in enumerate(result.tautomers, 1):
            is_canon = (t_smi == result.canonical_tautomer)
            canon_marker = "[bright_green]Yes[/bright_green]" if is_canon else "[dim]No[/dim]"
            smiles_style = "bright_cyan bold" if is_canon else "white"
            
            table.add_row(str(i), f"[{smiles_style}]{t_smi}[/{smiles_style}]", canon_marker)
            
        console.print(table)


# =============================================================================
# reaction command
# =============================================================================

@app.command(name="reaction")
def reaction_cmd(
    smiles: str = typer.Argument(..., help="Reaction SMILES (SMIRKS) string (e.g., A.B>>C)"),
) -> None:
    """Validate and analyze a Reaction SMILES (SMIRKS)."""
    result = validate_reaction(smiles)
    
    if not result.is_valid:
        console.print(f"[red]Error:[/red] {result.error}")
        raise typer.Exit(code=1)
        
    table = Table(title="Reaction SMILES Analysis", show_lines=True)
    table.add_column("Property", style="dim")
    table.add_column("Value", style="bright_cyan")
    
    table.add_row("Input", result.input_smiles)
    table.add_row("Reactants", str(result.num_reactants))
    table.add_row("Agents", str(result.num_agents))
    table.add_row("Products", str(result.num_products))
    
    console.print(table)
    console.print("[green]Reaction is perfectly valid![/green]")


# =============================================================================
# conformers command
# =============================================================================

@app.command(name="conformers")
def conformers_cmd(
    smiles: str = typer.Argument(..., help="Input SMILES string."),
    num: int = typer.Option(50, "--num-conformers", "--num", "-n", help="Number of conformers to generate."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output .sdf file (default: conformers.sdf)."),
) -> None:
    """Generate multiple 3D conformers for a molecule (ETKDG + MMFF)."""
    if output is None:
        output = Path("conformers.sdf")
    if output.suffix.lower() != ".sdf":
        console.print("[red]Error:[/red] Output file must have an .sdf extension.")
        raise typer.Exit(code=1)
        
    with console.status(f"[bold green]Generating {num} conformers...[/bold green]"):
        result = generate_conformers(smiles, num_conformers=num, output_sdf=str(output))
        
    if result.error:
        console.print(f"[red]Error:[/red] {result.error}")
        raise typer.Exit(code=1)
        
    console.print(f"\n[bright_green]Success![/bright_green] Generated {result.num_generated} optimized 3D conformers.")
    console.print(f"Saved to: [bright_cyan]{output}[/bright_cyan]")


# =============================================================================
# scaffold command
# =============================================================================

@app.command(name="scaffold")
def scaffold_cmd(
    smiles: str = typer.Argument(None, help="Input SMILES string."),
    file: Optional[Path] = typer.Option(None, "--file", "-f", help="Input file (CSV/SMI) for batch processing."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output CSV to save batch results."),
) -> None:
    """Extract the Murcko Scaffold framework from a SMILES string."""
    import csv
    if smiles is None and file is None:
        console.print("[red]Error:[/red] Provide either a SMILES argument or --file.")
        raise typer.Exit(code=1)

    if smiles and not file:
        result = extract_scaffold(smiles)
        if result.success:
            console.print(f"[dim]Input:[/dim]    {smiles}")
            if result.scaffold_smiles:
                console.print(f"[dim]Scaffold:[/dim] [bright_cyan]{escape(result.scaffold_smiles)}[/bright_cyan]")
            else:
                console.print(f"[dim]Scaffold:[/dim] [yellow]none — {result.note or 'acyclic molecule'}[/yellow]")
        else:
            console.print(f"[red]Error:[/red] {result.error}")
        return

    # Batch processing
    from chemlitmus.utils.parsers import parse_compounds_file
    compounds = parse_compounds_file(file)
    smiles_list = [c if isinstance(c, str) else c.get("smiles", "") for c in compounds]
    
    results = []
    with console.status(f"[bold green]Extracting scaffolds for {len(smiles_list)} compounds...[/bold green]"):
        for smi in smiles_list:
            if smi.strip():
                results.append(extract_scaffold(smi))
                
    table = Table(title=f"Scaffold Extraction ({len(results)} compounds)", show_lines=False)
    table.add_column("Input SMILES", style="dim", max_width=40)
    table.add_column("Murcko Scaffold", style="bright_cyan")
    
    for r in results:
        if not r.success:
            cell = f"[red]{escape(r.error or '')}[/red]"
        elif r.scaffold_smiles:
            cell = escape(r.scaffold_smiles)
        else:
            cell = f"[dim]— {escape(r.note or 'acyclic')}[/dim]"
        table.add_row(escape(r.input_smiles[:40]), cell)
    console.print(table)
    
    if output:
        with open(output, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["smiles", "scaffold", "error"])
            writer.writeheader()
            for r in results:
                writer.writerow({
                    "smiles": r.input_smiles,
                    "scaffold": r.scaffold_smiles or "",
                    "error": r.error or ""
                })
        console.print(f"[green]Saved scaffolds to:[/green] {output}")


# =============================================================================
# stereo command
# =============================================================================

@app.command(name="stereo")
def stereo_cmd(
    smiles: str = typer.Argument(..., help="Input SMILES string."),
    chiral_flag: bool = typer.Option(False, "--chiral-flag", help="Return error code 1 if stereochemistry is unassigned."),
) -> None:
    """Analyze stereocenters and unassigned stereochemistry."""
    result = analyze_stereochemistry(smiles)
    
    if not result.success:
        console.print(f"[red]Error:[/red] {result.error}")
        raise typer.Exit(code=1)
        
    console.print(f"\n[dim]Input SMILES:[/dim] {smiles}")
    
    if not result.chiral_centers:
        console.print("[dim]No stereocenters found.[/dim]")
        return
        
    table = Table(title="Stereocenters", show_lines=True)
    table.add_column("Atom Index", justify="right")
    table.add_column("Configuration", justify="center")
    
    for center in result.chiral_centers:
        conf = center["config"]
        color = "red" if conf == "?" else "bright_cyan"
        table.add_row(str(center["atom_idx"]), f"[{color}]{conf}[/{color}]")
        
    console.print(table)
    
    if result.has_unassigned:
        console.print("[yellow]Warning: Molecule has unassigned stereocenters ('?').[/yellow]")
        if chiral_flag:
            raise typer.Exit(code=1)
    else:
        console.print("[green]All stereocenters are fully assigned.[/green]")


# =============================================================================
# rgroup command
# =============================================================================

@app.command(name="rgroup")
def rgroup_cmd(
    core_arg: Optional[str] = typer.Argument(None, help="Core SMARTS (alternative to --core)."),
    core: Optional[str] = typer.Option(None, "--core", "-c", help="SMARTS string representing the core scaffold."),
    smiles: Optional[str] = typer.Option(None, "--smiles", "-s", help="Comma-separated SMILES strings to decompose."),
    file: Optional[Path] = typer.Option(None, "--file", "-f", help="CSV/SMI file containing SMILES."),
) -> None:
    """Perform R-Group Decomposition against a common core."""
    core = core or core_arg
    if not core:
        console.print("[red]Error:[/red] Provide the core SMARTS positionally or with --core.")
        raise typer.Exit(code=1)
    if not smiles and not file:
        console.print("[red]Error:[/red] Must provide either --smiles or --file.")
        raise typer.Exit(1)
        
    smiles_list = []
    if smiles:
        smiles_list.extend([s.strip() for s in smiles.split(",") if s.strip()])
        
    if file:
        from chemlitmus.utils.parsers import parse_compounds_file
        compounds = parse_compounds_file(file)
        smiles_list.extend([c if isinstance(c, str) else c.get("smiles", "") for c in compounds])
        
    with console.status(f"[bold green]Decomposing {len(smiles_list)} molecules...[/bold green]"):
        results = rgroup_decomposition(core, smiles_list)
        
    table = Table(title="R-Group Decomposition", show_lines=True)
    table.add_column("Input SMILES", style="dim")
    
    # Collect all unique R-group labels
    all_keys = set()
    for r in results:
        if r.is_matched:
            all_keys.update(r.decomposition.keys())
            
    # Sort keys: Core first, then R1, R2, etc.
    sorted_keys = sorted(list(all_keys), key=lambda k: (0 if k == "Core" else 1, k))
    
    for k in sorted_keys:
        table.add_column(k, justify="center")
        
    table.add_column("Status", justify="right")
    
    for r in results:
        row = [r.input_smiles]
        for k in sorted_keys:
            if r.is_matched:
                row.append(r.decomposition.get(k, ""))
            else:
                row.append("")
        if r.is_matched:
            row.append("[green]Matched[/green]")
        else:
            row.append(f"[red]{r.error}[/red]")
        table.add_row(*row)
        
    console.print(table)


# =============================================================================
# augment command
# =============================================================================

@app.command(name="augment")
def augment_cmd(
    smiles: str = typer.Argument(..., help="Input SMILES string to augment."),
    num: int = typer.Option(5, "--num", "-n", help="Number of augmented SMILES to generate."),
) -> None:
    """Generate uncanonical/randomized SMILES strings for data augmentation."""
    result = augment_smiles(smiles, num_augmentations=num)
    
    if not result.success:
        console.print(f"[red]Error:[/red] {result.error}")
        raise typer.Exit(1)
        
    console.print(f"\n[dim]Input:[/dim] {smiles}")
    console.print(f"[bold green]Generated {len(result.augmented_smiles)} variations:[/bold green]")
    
    for i, s in enumerate(result.augmented_smiles, 1):
        console.print(f"  {i}. [bright_cyan]{s}[/bright_cyan]")


# =============================================================================
# atommap command
# =============================================================================

@app.command(name="atommap")
def atommap_cmd(
    smiles: str = typer.Argument(..., help="Input SMILES string to map."),
) -> None:
    """Assign unique atom map numbers to all atoms in a molecule."""
    result = map_atoms(smiles)
    
    if not result.success:
        console.print(f"[red]Error:[/red] {result.error}")
        raise typer.Exit(1)
        
    console.print(f"\n[dim]Input:[/dim]  {smiles}")
    console.print(f"[dim]Mapped:[/dim] [bright_cyan]{result.mapped_smiles}[/bright_cyan]\n")

# =============================================================================
# validate command
# =============================================================================

@app.command(name="validate")
def validate_cmd(
    smiles: str = typer.Argument(None, help="SMILES string to validate."),
    file: Optional[Path] = typer.Option(None, "--file", "-f", help="Batch input (CSV/TSV/XLSX/SMI/SDF)."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Save batch results to CSV."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Print only the canonical SMILES (single mode) - useful in pipelines."),
) -> None:
    """Validate SMILES offline and report canonical form plus basic properties.

    Exit code 0 if valid, 2 if invalid (single mode), so it works as a shell guard.
    For an explanation of *why* a SMILES fails, use `chemlitmus diagnose`.
    """
    import csv

    if smiles is None and file is None:
        console.print("[red]Error:[/red] Provide a SMILES argument or --file.")
        raise typer.Exit(code=1)

    if smiles is not None and file is None:
        r = validate_smiles(smiles)
        if not r.is_valid:
            if not quiet:
                console.print(f"[red]Invalid:[/red] {r.error_message}")
                console.print("[dim]Run `chemlitmus diagnose` on this string for a located explanation.[/dim]")
            raise typer.Exit(code=2)
        if quiet:
            print(r.canonical_smiles)
            return
        tbl = Table(title="[bold]SMILES Validation[/bold]", show_header=True, header_style="bold magenta")
        tbl.add_column("Property", style="cyan"); tbl.add_column("Value")
        tbl.add_row("Input", escape(r.input_smiles))
        tbl.add_row("Canonical SMILES", escape(r.canonical_smiles or ""))
        tbl.add_row("Formula", r.molecular_formula or "")
        tbl.add_row("Exact mass", f"{r.molecular_weight:.4f}" if r.molecular_weight is not None else "")
        tbl.add_row("LogP", f"{r.logp:.2f}" if r.logp is not None else "")
        tbl.add_row("H-bond donors / acceptors", f"{r.hbd} / {r.hba}")
        tbl.add_row("TPSA", f"{r.tpsa:.2f}" if r.tpsa is not None else "")
        tbl.add_row("Heavy atoms", str(r.heavy_atom_count))
        console.print(tbl)
        return

    records = parse_compounds_file(file)
    results = [validate_smiles(s) for s in records]
    n_ok = sum(1 for r in results if r.is_valid)
    tbl = Table(title=f"[bold]SMILES Validation — {len(results)} records[/bold]", show_header=True, header_style="bold magenta")
    tbl.add_column("Input", max_width=36, overflow="fold"); tbl.add_column("Valid", justify="center")
    tbl.add_column("Canonical", max_width=36, overflow="fold"); tbl.add_column("MW", justify="right")
    for r in results[:20]:
        tbl.add_row(escape(r.input_smiles), "[green]yes[/green]" if r.is_valid else "[red]no[/red]",
                    escape(r.canonical_smiles or ""), f"{r.molecular_weight:.2f}" if r.molecular_weight is not None else "")
    console.print(tbl)
    if len(results) > 20:
        console.print(f"[dim]... {len(results) - 20} more; use --output to save all.[/dim]")
    console.print(f"\n[bold]Summary[/bold]: valid={n_ok}  invalid={len(results) - n_ok}")
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["input_smiles", "is_valid", "canonical_smiles", "molecular_formula",
                                               "molecular_weight", "logp", "hbd", "hba", "tpsa", "heavy_atom_count", "error_message"])
            w.writeheader()
            for r in results:
                w.writerow({k: ("" if v is None else v) for k, v in r.model_dump().items()})
        console.print(f"[green]Saved:[/green] {output}")


# =============================================================================
# identity command
# =============================================================================

@app.command(name="identity")
def identity_cmd(
    smiles: str = typer.Argument(None, help="Single SMILES: show its identity keys at every level."),
    file: Optional[Path] = typer.Option(None, "--file", "-f", help="Collection to group (CSV/.smi/.sdf)."),
    level: str = typer.Option("parent", "--level", "-l", help="Identity level: exact, parent, tautomer, nostereo, skeleton, formula."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Write per-record identity keys and group ids to CSV."),
    show: int = typer.Option(15, "--show", help="Number of multi-member groups to list."),
) -> None:
    """Compute layered molecular identity, or group a collection by identity level.

    Levels nest from strictest to loosest: exact (salts and charges included) > parent
    (largest fragment, neutralised) > tautomer / nostereo > skeleton > formula. Grouping a
    collection at a level shows which records are the same compound at that resolution and
    what varies among them (salt form, tautomer, stereochemistry).
    """
    import csv

    if level not in IDENTITY_LEVELS:
        console.print(f"[red]Error:[/red] Unknown level {level!r}. Valid: {', '.join(IDENTITY_LEVELS)}")
        raise typer.Exit(code=1)
    if smiles is None and file is None:
        console.print("[red]Error:[/red] Provide a SMILES argument or --file.")
        raise typer.Exit(code=1)

    if smiles is not None and file is None:
        k = compute_identity(smiles)
        if not k.is_valid:
            console.print(f"[red]Error:[/red] {k.error}")
            raise typer.Exit(code=1)
        t = Table(title="[bold]Identity Keys[/bold]", show_header=True, header_style="bold magenta")
        t.add_column("Level", style="cyan"); t.add_column("Key")
        for lv in IDENTITY_LEVELS:
            t.add_row(lv, escape(k.key(lv) or ""))
        console.print(t)
        console.print(f"[dim]fragments={k.n_fragments}  charged={'yes' if k.had_charge else 'no'}  stereo specified={'yes' if k.has_stereo else 'no'}[/dim]")
        return

    records = parse_compounds_file(file)
    with console.status(f"[bold green]Computing identity for {len(records)} records...[/bold green]"):
        rep = group_by_identity(records, level=level)

    t = Table(title=f"[bold]Identity Report — {rep.n_records} records at level '{level}'[/bold]", show_header=True, header_style="bold magenta")
    t.add_column("Level", style="cyan"); t.add_column("Distinct compounds", justify="right"); t.add_column("Collapsed records", justify="right")
    for lv in IDENTITY_LEVELS:
        n = rep.n_groups_by_level.get(lv, 0)
        mark = " [bold]<[/bold]" if lv == level else ""
        t.add_row(lv + mark, str(n), str(rep.n_valid - n))
    console.print(t)
    console.print(f"[dim]valid={rep.n_valid}  invalid={rep.n_records - rep.n_valid}  distinct at '{level}'={rep.n_groups}  removable duplicates={rep.n_collapsed}[/dim]")

    if rep.groups:
        g = Table(title=f"Groups with >1 member (first {min(show, len(rep.groups))} of {len(rep.groups)})", show_header=True, header_style="bold blue")
        g.add_column("Size", justify="right"); g.add_column("Varies by", style="yellow"); g.add_column("Members", max_width=70, overflow="fold")
        for grp in rep.groups[:show]:
            g.add_row(str(grp.size), ", ".join(grp.differs_by) or "representation only", escape("  |  ".join(grp.smiles[:4]) + ("  ..." if grp.size > 4 else "")))
        console.print(g)

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        key_to_group = {}
        for grp in rep.groups:
            key_to_group[grp.key] = grp
        with open(output, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["index", "input_smiles", "is_valid"] + IDENTITY_LEVELS + [f"group_id_{level}", f"group_size_{level}", "group_varies_by"])
            group_ids = {}
            for i, k in enumerate(rep.keys):
                key = k.key(level) if k.is_valid else None
                gid = group_ids.setdefault(key, len(group_ids)) if key is not None else ""
                grp = key_to_group.get(key)
                w.writerow([i, k.input_smiles, k.is_valid] + [(k.key(lv) or "") for lv in IDENTITY_LEVELS]
                           + [gid, grp.size if grp else (1 if key else ""), ";".join(grp.differs_by) if grp else ""])
        console.print(f"[green]Saved:[/green] {output}")


# =============================================================================
# diff command
# =============================================================================

@app.command(name="diff")
def diff_cmd(
    file_a: Path = typer.Argument(..., help="Reference (older) collection."),
    file_b: Path = typer.Argument(..., help="Comparison (newer) collection."),
    level: str = typer.Option("parent", "--level", "-l", help="Identity level for matching: exact, parent, tautomer, nostereo, skeleton, formula."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Write added/removed/changed entries to CSV."),
    json_out: Optional[Path] = typer.Option(None, "--json", help="Write the complete result to JSON."),
    include_unchanged: bool = typer.Option(False, "--include-unchanged", help="Also list unchanged compounds."),
    show: int = typer.Option(10, "--show", help="Entries to list per category in the terminal."),
) -> None:
    """Structure-aware comparison of two compound collections.

    Matches compounds at the chosen identity level and reports what was added, removed, left
    unchanged, or kept but written differently (salt form, tautomer, stereochemistry), so
    database releases, vendor catalogues and generated libraries can be compared without
    being fooled by SMILES formatting.
    """
    import csv
    import json

    if level not in IDENTITY_LEVELS:
        console.print(f"[red]Error:[/red] Unknown level {level!r}. Valid: {', '.join(IDENTITY_LEVELS)}")
        raise typer.Exit(code=1)
    a, b = parse_compounds_file(file_a), parse_compounds_file(file_b)
    with console.status(f"[bold green]Comparing {len(a)} vs {len(b)} records at level '{level}'...[/bold green]"):
        res = diff_libraries(a, b, level=level, include_unchanged=include_unchanged)

    t = Table(title=f"[bold]Library Diff — {file_a.name} → {file_b.name} at level '{level}'[/bold]", show_header=True, header_style="bold magenta")
    t.add_column("", style="cyan"); t.add_column("Count", justify="right"); t.add_column("Meaning", style="dim")
    t.add_row("Compounds in A", str(res.n_keys_a), f"{res.n_valid_a} valid records")
    t.add_row("Compounds in B", str(res.n_keys_b), f"{res.n_valid_b} valid records")
    t.add_row("[green]Added[/green]", str(res.n_added), "in B only")
    t.add_row("[red]Removed[/red]", str(res.n_removed), "in A only")
    t.add_row("Unchanged", str(res.n_unchanged), "same compound, same representation")
    t.add_row("[yellow]Changed[/yellow]", str(res.n_changed), "same compound at this level, different representation")
    t.add_row("Multiplicity changes", str(res.multiplicity_changes), "record count differs between A and B")
    t.add_row("Overlap (Jaccard)", f"{res.jaccard:.3f}", "shared compounds / all compounds")
    console.print(t)
    if res.changes_by_kind:
        k = Table(title="What changed", show_header=True, header_style="bold blue")
        k.add_column("Kind", style="yellow"); k.add_column("Compounds", justify="right")
        for kind, n in sorted(res.changes_by_kind.items(), key=lambda x: -x[1]):
            k.add_row(kind, str(n))
        console.print(k)

    for status, style in (("removed", "red"), ("added", "green"), ("changed", "yellow")):
        rows = [e for e in res.entries if e.status == status]
        if not rows:
            continue
        et = Table(title=f"[{style}]{status.capitalize()}[/{style}] (first {min(show, len(rows))} of {len(rows)})", show_header=True, header_style="bold blue")
        et.add_column("A", max_width=40, overflow="fold"); et.add_column("B", max_width=40, overflow="fold")
        if status == "changed":
            et.add_column("Change", style="yellow")
        for e in rows[:show]:
            cells = [escape(" | ".join(e.smiles_a)), escape(" | ".join(e.smiles_b))]
            if status == "changed":
                cells.append(e.change or "")
            et.add_row(*cells)
        console.print(et)

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["status", "level", "key", "smiles_a", "smiles_b", "count_a", "count_b", "change"])
            w.writeheader()
            for e in res.entries:
                w.writerow({"status": e.status, "level": level, "key": e.key, "smiles_a": "|".join(e.smiles_a),
                            "smiles_b": "|".join(e.smiles_b), "count_a": e.count_a, "count_b": e.count_b, "change": e.change or ""})
        console.print(f"[green]Saved:[/green] {output}")
    if json_out:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(json.dumps(res.model_dump(), indent=2))
        console.print(f"[green]Saved:[/green] {json_out}")


# =============================================================================
# diagnose command
# =============================================================================

def _print_diagnosis(d: SmilesDiagnosis) -> None:
    if d.is_valid:
        console.print(Panel(Text(f"  {d.input_smiles}\n  valid  →  {d.canonical_smiles}", style="green"), title="[bold]SMILES Diagnosis[/bold]", border_style="green"))
        return
    body = Text()
    body.append(f"  {d.input_smiles}\n", style="bold")
    body.append(f"  {d.caret_line()}\n", style="bold red")
    for p in d.problems:
        where = f"position {p.position}" if p.position is not None else ("atom " + str(p.atom_index) if p.atom_index is not None else "")
        body.append(f"  [{p.category}] ", style="yellow"); body.append(f"{p.message}" + (f"  ({where})" if where else "") + "\n")
        if p.suggestion:
            body.append(f"      → {p.suggestion}\n", style="dim")
    if d.repaired_smiles is not None:
        if d.repaired_is_valid:
            body.append(f"\n  Repaired: {d.repaired_smiles}\n", style="green")
        else:
            body.append(f"\n  Attempted repair still invalid: {d.repaired_smiles}\n", style="red")
        for r in d.repairs_applied:
            body.append(f"      • {r}\n", style="dim")
    console.print(Panel(body, title="[bold]SMILES Diagnosis[/bold]", border_style="red"))


@app.command(name="diagnose")
def diagnose_cmd(
    smiles: str = typer.Argument(None, help="SMILES string to diagnose."),
    file: Optional[Path] = typer.Option(None, "--file", "-f", help="Batch input (CSV/.smi)."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Write per-record diagnoses to CSV."),
    no_repair: bool = typer.Option(False, "--no-repair", help="Do not attempt mechanical repairs."),
    only_invalid: bool = typer.Option(True, "--only-invalid/--all", help="In batch mode, report only invalid records (default) or all."),
) -> None:
    """Explain why a SMILES string fails to parse, with character positions and suggested fixes.

    Checks run in a fixed order — characters, bracket atoms, parentheses, ring closures,
    RDKit syntax, valence, aromaticity — and each finding points at a position or an atom.
    Safe mechanical repairs (whitespace, dangling branches, unclosed ring digits, [nH],
    non-ring aromatic atoms) are attempted and re-validated.
    """
    import csv

    if smiles is None and file is None:
        console.print("[red]Error:[/red] Provide a SMILES argument or --file.")
        raise typer.Exit(code=1)

    if smiles is not None and file is None:
        d = diagnose_smiles(smiles, try_repair=not no_repair)
        _print_diagnosis(d)
        raise typer.Exit(code=0 if d.is_valid else 2)

    records = parse_compounds_file(file)
    results = []
    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"), BarColumn(), TaskProgressColumn(), console=console) as progress:
        task = progress.add_task("Diagnosing...", total=len(records))
        for s in records:
            results.append(diagnose_smiles(s, try_repair=not no_repair))
            progress.advance(task)
    invalid = [d for d in results if not d.is_valid]
    repaired = sum(1 for d in invalid if d.repaired_is_valid)
    by_cat = {}
    for d in invalid:
        c = d.primary_category or "unknown"
        by_cat[c] = by_cat.get(c, 0) + 1

    t = Table(title=f"[bold]SMILES Diagnosis — {len(records)} records[/bold]", show_header=True, header_style="bold magenta")
    t.add_column("Primary problem", style="cyan"); t.add_column("Records", justify="right")
    for c, n in sorted(by_cat.items(), key=lambda x: -x[1]):
        t.add_row(c, str(n))
    t.add_row("[bold]invalid total[/bold]", f"[bold]{len(invalid)}[/bold]")
    t.add_row("repaired mechanically", str(repaired))
    t.add_row("valid", str(len(records) - len(invalid)))
    console.print(t)

    shown = invalid if only_invalid else results
    for d in shown[:10]:
        _print_diagnosis(d)
    if len(shown) > 10:
        console.print(f"[dim]... {len(shown) - 10} more; use --output to save all.[/dim]")

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["input_smiles", "is_valid", "canonical_smiles", "primary_category", "n_problems",
                                               "problems", "positions", "suggestions", "repaired_smiles", "repaired_is_valid", "repairs_applied"])
            w.writeheader()
            for d in (results if not only_invalid else invalid):
                w.writerow({
                    "input_smiles": d.input_smiles, "is_valid": d.is_valid, "canonical_smiles": d.canonical_smiles or "",
                    "primary_category": d.primary_category or "", "n_problems": len(d.problems),
                    "problems": " | ".join(f"[{p.category}] {p.message}" for p in d.problems),
                    "positions": ";".join("" if p.position is None else str(p.position) for p in d.problems),
                    "suggestions": " | ".join(p.suggestion for p in d.problems if p.suggestion),
                    "repaired_smiles": d.repaired_smiles or "", "repaired_is_valid": "" if d.repaired_is_valid is None else d.repaired_is_valid,
                    "repairs_applied": " | ".join(d.repairs_applied),
                })
        console.print(f"[green]Saved:[/green] {output}")


# =============================================================================
# smartsaudit command
# =============================================================================

def _print_audit_summary(res: SmartsAuditResult, breadth_threshold: float) -> None:
    n = res.n_patterns
    def pct(k: int) -> str:
        return f"{100.0 * k / n:.1f}%" if n else "-"

    t = Table(title=f"[bold]SMARTS Audit — {n} patterns × {res.n_molecules:,} reference molecules[/bold]",
              show_header=True, header_style="bold magenta")
    t.add_column("Check", style="cyan")
    t.add_column("Patterns", justify="right")
    t.add_column("Share", justify="right")
    t.add_column("Meaning", style="dim")
    t.add_row("Unparseable", str(res.n_unparseable), pct(res.n_unparseable), "RDKit cannot compile the SMARTS")
    t.add_row("Need explicit H", str(res.n_needs_explicit_h), pct(res.n_needs_explicit_h),
              "Contain a hydrogen atom; silently dead under the default preparation")
    if "breadth" in res.checks_run:
        t.add_row(f"Over-broad (>{breadth_threshold:.0%})", str(res.n_over_broad), pct(res.n_over_broad),
                  "Match a large share of the reference set")
    if "dead" in res.checks_run:
        t.add_row("Never fire", str(res.n_dead), pct(res.n_dead), "Zero hits on the reference set")
        t.add_row("  ↳ never-matching atom", str(res.n_dead_never_matching_atom), pct(res.n_dead_never_matching_atom),
                  "A query atom matches no real atom — likely defective")
    if "redundancy" in res.checks_run:
        t.add_row("Exact duplicates", str(res.n_duplicates), pct(res.n_duplicates), "Identical SMARTS string appears earlier")
        t.add_row("Library-equivalent", str(res.n_equivalent), pct(res.n_equivalent), "Identical hit set to another pattern")
        t.add_row("Strictly subsumed", str(res.n_subsumed), pct(res.n_subsumed), "Another pattern's hits contain all of these")
    if res.sensitivity is not None:
        t.add_row("Preparation-sensitive", str(res.sensitivity.n_sensitive_patterns), pct(res.sensitivity.n_sensitive_patterns),
                  "Hit count changes with molecule preparation")
    t.add_row("[bold]No flags[/bold]", f"[bold]{res.n_clean}[/bold]", f"[bold]{pct(res.n_clean)}[/bold]", "")
    console.print(t)

    if res.sensitivity is not None:
        s = res.sensitivity
        st = Table(title="[bold]Reproducibility across molecule preparations[/bold]", show_header=True, header_style="bold blue")
        st.add_column("Preparation", style="cyan")
        st.add_column("Compounds flagged", justify="right")
        st.add_column("Share", justify="right")
        st.add_column("Patterns firing", justify="right")
        st.add_column("Verdict flips vs default", justify="right")
        for prep in s.preparations:
            flagged = s.compounds_flagged.get(prep, 0)
            flips = s.verdict_flips.get(prep)
            st.add_row(prep, f"{flagged:,}", f"{100.0*flagged/s.n_molecules:.1f}%" if s.n_molecules else "-",
                       str(s.patterns_firing.get(prep, 0)),
                       "—" if flips is None else f"{flips:,} ({100.0*flips/s.n_molecules:.1f}%)")
        console.print(st)


def _print_explanation(e: SmartsExplanation) -> None:
    body = Text()
    body.append(f"  SMARTS:      {e.smarts}\n", style="cyan")  # Text(): no markup parsing
    if not e.parses:
        body.append(f"  Error:       {e.parse_error}\n", style="red")
        console.print(Panel(body, title="[bold]SMARTS Explanation[/bold]", border_style="red"))
        return
    body.append(f"  Normalized:  {e.normalized_smarts}\n", style="dim")
    body.append(f"  Atoms/bonds: {e.n_query_atoms} / {e.n_query_bonds}\n")
    body.append(f"  Recursive:   {'yes' if e.has_recursive_smarts else 'no'}      Needs explicit H: {'yes' if e.requires_explicit_h else 'no'}\n")
    body.append("  Hits:        " + "   ".join(f"{k}={v:,}" for k, v in e.hits_by_preparation.items()) + f"   (of {e.n_molecules:,})\n")
    verdict_style = "green" if e.verdict == "fires" else ("yellow" if e.verdict.startswith("fires only") else "red")
    body.append(f"  Verdict:     {e.verdict}\n", style=f"bold {verdict_style}")
    console.print(Panel(body, title="[bold]SMARTS Explanation[/bold]", border_style="cyan"))

    at = Table(title="Atom primitives (each tested alone against a reference sample)", show_header=True, header_style="bold blue")
    at.add_column("#", justify="right")
    at.add_column("Primitive", style="cyan")
    at.add_column("Molecules with a matching atom", justify="right")
    for a in e.atoms:
        n = f"[red]{a.n_matching_molecules}[/red]" if a.n_matching_molecules == 0 else str(a.n_matching_molecules)
        at.add_row(str(a.atom_index), escape(a.query), n)
    console.print(at)
    if e.example_matches:
        console.print("[bold]Example matches:[/bold]")
        for s in e.example_matches:
            console.print(f"  [green]{escape(s)}[/green]")


@app.command(name="smartsaudit")
def smartsaudit_cmd(
    patterns: Optional[Path] = typer.Argument(None, help="Pattern file: CSV/TSV/XLSX with a 'smarts' column, or a text file with one SMARTS per line."),
    explain: Optional[str] = typer.Option(None, "--explain", "-e", help="Explain a single SMARTS instead of auditing a file."),
    library: Optional[Path] = typer.Option(None, "--library", "-l", help="Reference molecules (.smi/.csv/.sdf). Default: bundled ChEMBL-derived set."),
    checks: str = typer.Option("all", "--checks", "-c", help="Comma-separated subset of: compile,breadth,dead,redundancy,sensitivity (or all)."),
    breadth_threshold: float = typer.Option(0.10, "--breadth-threshold", help="Flag patterns matching more than this fraction of the reference set."),
    max_molecules: Optional[int] = typer.Option(None, "--max-molecules", help="Use only the first N reference molecules."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Write the per-pattern audit table to CSV."),
    json_out: Optional[Path] = typer.Option(None, "--json", help="Write the complete result (including the sensitivity summary) to JSON."),
    show: int = typer.Option(15, "--show", help="Number of flagged patterns to list in the terminal."),
) -> None:
    """Audit a SMARTS pattern set (structural alerts, substructure filters) for defects.

    Reports patterns that fail to compile, need explicit hydrogens, match too broadly, never
    match anything (triaged into never-matching-atom vs rare-combination), duplicate or subsume
    one another, and change their hit count with molecule preparation (implicit H, explicit H,
    kekulized). Fully offline; matching runs multithreaded via RDKit's SubstructLibrary.
    """
    import csv
    import json

    if explain is None and patterns is None:
        console.print("[red]Error:[/red] Provide a pattern file, or --explain SMARTS.")
        raise typer.Exit(code=1)

    with console.status("[bold green]Loading reference library...[/bold green]"):
        try:
            mols, source = load_reference_library(library, max_molecules=max_molecules)
        except Exception as exc:
            console.print(f"[red]Error loading reference library:[/red] {exc}")
            raise typer.Exit(code=1)

    if explain is not None:
        with console.status("[bold green]Explaining pattern...[/bold green]"):
            e = explain_smarts(explain, library=mols)
        _print_explanation(e)
        if json_out:
            json_out.parent.mkdir(parents=True, exist_ok=True)
            json_out.write_text(json.dumps(e.model_dump(), indent=2))
            console.print(f"[green]Saved:[/green] {json_out}")
        return

    try:
        triples = load_patterns(patterns)
    except Exception as exc:
        console.print(f"[red]Error reading patterns:[/red] {exc}")
        raise typer.Exit(code=1)
    if not triples:
        console.print("[red]Error:[/red] No patterns found in file.")
        raise typer.Exit(code=1)

    check_list = None if checks.strip().lower() == "all" else [c for c in checks.split(",") if c.strip()]
    with console.status(f"[bold green]Auditing {len(triples)} patterns against {len(mols):,} molecules...[/bold green]"):
        try:
            res = audit_smarts(triples, library=mols, library_source=source, checks=check_list,
                               breadth_threshold=breadth_threshold)
        except ValueError as exc:
            console.print(f"[red]Error:[/red] {exc}")
            raise typer.Exit(code=1)

    console.print(f"[dim]Reference library: {source}  ·  {res.elapsed_seconds:.1f}s[/dim]")
    _print_audit_summary(res, breadth_threshold)

    flagged = [p for p in res.patterns if not p.clean]
    if flagged:
        ft = Table(title=f"Flagged patterns (first {min(show, len(flagged))} of {len(flagged)})", show_header=True, header_style="bold magenta")
        ft.add_column("#", justify="right")
        ft.add_column("Name", style="cyan", max_width=22, overflow="fold")
        ft.add_column("SMARTS", max_width=44, overflow="fold")
        ft.add_column("Hits", justify="right")
        ft.add_column("Flags", style="yellow")
        for p in flagged[:show]:
            ft.add_row(str(p.index), escape(p.name or ""), escape(p.smarts), "" if p.n_hits is None else str(p.n_hits), ", ".join(p.flags))
        console.print(ft)

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        rows = res.to_rows()
        with open(output, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        console.print(f"[green]Per-pattern table saved to:[/green] {output}")
    if json_out:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(json.dumps(res.model_dump(), indent=2))
        console.print(f"[green]Full result saved to:[/green] {json_out}")


# =============================================================================
# update command
# =============================================================================

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


if __name__ == "__main__":
    app()
