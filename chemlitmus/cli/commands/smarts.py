"""SMARTS catalogue quality control commands."""

from pathlib import Path
from typing import Optional
import typer
from rich.table import Table
from rich.panel import Panel
from rich.text import Text
from rich.markup import escape

from chemlitmus import (
    audit_smarts, explain_smarts, load_patterns, load_reference_library,
    SmartsAuditResult, SmartsExplanation,
)
from chemlitmus.cli._app import app, console


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
