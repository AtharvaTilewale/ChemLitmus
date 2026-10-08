"""Standardisation, naming and molecular identity commands."""

from pathlib import Path
from typing import Optional
import typer
from rich.table import Table
from rich.markup import escape

from chemlitmus import (
    standardize_smiles, StandardizeResult, STANDARDIZE_STEPS,
    get_iupac_name, IUPACResult,
    enumerate_tautomers, TautomerResult,
    analyze_stereochemistry, compute_identity, group_by_identity, IDENTITY_LEVELS,
    diff_libraries,
)
from chemlitmus.utils.parsers import parse_compounds_file
from chemlitmus.cli._app import app, console


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
