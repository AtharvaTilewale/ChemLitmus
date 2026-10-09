"""Validation and diagnosis commands."""

from pathlib import Path
from typing import Optional
import typer
from rich.table import Table
from rich.panel import Panel
from rich.text import Text
from rich.markup import escape
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from chemlitmus import (
    validate_smiles,
    diagnose_smiles, SmilesDiagnosis,
)
from chemlitmus.utils.parsers import parse_compounds_file
from chemlitmus.cli._app import app, console


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
            body.append(f"\n  Candidate repair (parses; not verified to be the intended molecule): {d.repaired_smiles}\n", style="green")
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
    t.add_row("candidate repairs (mechanical, unverified)", str(repaired))
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
                                               "problems", "positions", "suggestions", "repaired_smiles", "repaired_is_valid", "repair_status", "repairs_applied"])
            w.writeheader()
            for d in (results if not only_invalid else invalid):
                w.writerow({
                    "input_smiles": d.input_smiles, "is_valid": d.is_valid, "canonical_smiles": d.canonical_smiles or "",
                    "primary_category": d.primary_category or "", "n_problems": len(d.problems),
                    "problems": " | ".join(f"[{p.category}] {p.message}" for p in d.problems),
                    "positions": ";".join("" if p.position is None else str(p.position) for p in d.problems),
                    "suggestions": " | ".join(p.suggestion for p in d.problems if p.suggestion),
                    "repaired_smiles": d.repaired_smiles or "", "repaired_is_valid": "" if d.repaired_is_valid is None else d.repaired_is_valid, "repair_status": d.repair_status,
                    "repairs_applied": " | ".join(d.repairs_applied),
                })
        console.print(f"[green]Saved:[/green] {output}")
