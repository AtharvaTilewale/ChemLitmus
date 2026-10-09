"""SMARTS catalogue quality control commands."""

from pathlib import Path
from typing import List, Optional, Tuple
import typer
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn
from rich.panel import Panel
from rich.text import Text
from rich.markup import escape

from chemlitmus import (
    audit_smarts, explain_smarts, load_patterns, load_reference_library,
    SmartsAuditResult, SmartsExplanation,
)
from chemlitmus.cli._app import app, console
from rdkit import Chem
from chemlitmus.core.smartsaudit import AUDIT_CHECKS


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
    if res.reference_panel:
        p = res.reference_panel
        console.print(f"  [dim]Reference panel: {p.n_molecules} molecules ({escape(p.source)}), median {p.median_heavy_atoms:g} heavy atoms, "
                      f"{p.fraction_with_ring:.0%} with a ring, {p.fraction_charged:.0%} charged; top elements " + ", ".join(f"{k} {v}" for k, v in list(p.elements.items())[:6]) + "[/dim]")
    if res.catalogue and (res.catalogue.name or len(res.catalogue.rule_sets) > 1):
        c = res.catalogue
        console.print(f"  [dim]Catalogue: {escape(c.name or 'unnamed')}" + (f" v{escape(c.version)}" if c.version else "") + (f" · {escape(c.licence)}" if c.licence else "")
                      + f" · {c.n_patterns} patterns in {len(c.rule_sets)} set(s): " + ", ".join(f"{escape(k)} {v}" for k, v in list(c.rule_sets.items())[:6])
                      + " — measurements cover every set listed, not one named set[/dim]")
    for h in res.holdout:
        console.print(f"  [yellow]Holdout {escape(h.panel)}:[/yellow] {h.n_patterns_firing} pattern(s) fire on {h.n_molecules} molecules; "
                      f"{len(h.patterns_revived)} of them never fired on the reference panel — 'not observed' is panel-specific.")
    if res.match_semantics:
        ms = res.match_semantics
        console.print(f"  [dim]Match semantics: preparation {ms.preparation}, hydrogens {ms.hydrogens}, chirality {'on' if ms.use_chirality else 'ignored (RDKit default)'}, RDKit {escape(ms.rdkit_version or '?')}[/dim]")
    if "proof" in res.checks_run:
        t.add_row("Proven redundant", str(res.n_proven_redundant), pct(res.n_proven_redundant), "Static proof: another pattern contains it, for every molecule")
        t.add_row("Proof undecided", str(res.n_proof_undecided), pct(res.n_proof_undecided), "No witness found — not refuted")
        t.add_row("Not provable", str(res.n_proof_not_analysable), pct(res.n_proof_not_analysable), "Recursive SMARTS, isotope, valence, ring size, disconnected")
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
    checks: str = typer.Option("all", "--checks", "-c", help="Comma-separated subset of: compile,breadth,dead,redundancy,sensitivity,proof (or all; proof is opt-in, e.g. all,proof)."),
    breadth_threshold: float = typer.Option(0.10, "--breadth-threshold", help="Flag patterns matching more than this fraction of the reference set."),
    max_molecules: Optional[int] = typer.Option(None, "--max-molecules", help="Use only the first N reference molecules."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Write the per-pattern audit table to CSV."),
    json_out: Optional[Path] = typer.Option(None, "--json", help="Write the complete result (including the sensitivity summary) to JSON."),
    show: int = typer.Option(15, "--show", help="Number of flagged patterns to list in the terminal."),
    holdout: Optional[List[Path]] = typer.Option(None, "--holdout", help="Additional library to test patterns against; a 'never fires' verdict is panel-specific. Repeatable."),
    catalogue_name: Optional[str] = typer.Option(None, "--catalogue-name", help="Catalogue provenance recorded in the result."),
    catalogue_version: Optional[str] = typer.Option(None, "--catalogue-version"),
    catalogue_source: Optional[str] = typer.Option(None, "--catalogue-source", help="URL or citation of the catalogue as distributed."),
    catalogue_licence: Optional[str] = typer.Option(None, "--catalogue-licence"),
    cleanup: Optional[Path] = typer.Option(None, "--cleanup", help="Write reviewable cleanup proposals to this CSV. Nothing is ever removed automatically."),
    allow_cross_set: bool = typer.Option(False, "--allow-cross-set-cleanup", help="Also recommend removing a rule whose only cover is in a different published set (loses provenance)."),
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

    check_list = None if checks.strip().lower() == "all" else [c.strip().lower() for c in checks.split(",") if c.strip()]
    from chemlitmus.core.smartsaudit import CatalogueMetadata
    cat = CatalogueMetadata(name=catalogue_name, version=catalogue_version, source=catalogue_source, licence=catalogue_licence) \
        if any((catalogue_name, catalogue_version, catalogue_source, catalogue_licence)) else None
    holdout_libs = {}
    for h in holdout or []:
        if not h.exists():
            console.print(f"[red]Error:[/red] Holdout library not found: {h}")
            raise typer.Exit(code=1)
        hm, _src = load_reference_library(h)
        holdout_libs[h.stem] = hm
    if check_list and any(c not in AUDIT_CHECKS and c != "all" for c in check_list):
        bad = [c for c in check_list if c not in AUDIT_CHECKS and c != "all"]
        console.print(f"[red]Error:[/red] Unknown check(s) {bad}. Valid: {', '.join(AUDIT_CHECKS)} or all.")
        raise typer.Exit(code=1)
    with console.status(f"[bold green]Auditing {len(triples)} patterns against {len(mols):,} molecules...[/bold green]"):
        try:
            res = audit_smarts(triples, library=mols, library_source=source, checks=check_list,
                               breadth_threshold=breadth_threshold,
                               catalogue=cat, holdout_libraries=holdout_libs or None)
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
    if cleanup:
        import csv as _csv

        from chemlitmus.core.cleanup import propose_cleanup

        crep = propose_cleanup(res, allow_cross_set=allow_cross_set)
        cleanup.parent.mkdir(parents=True, exist_ok=True)
        with open(cleanup, "w", newline="", encoding="utf-8") as fh:
            w = _csv.writer(fh)
            w.writerow(["index", "name", "rule_set", "smarts", "action", "evidence", "covered_by", "covered_by_rule_sets", "crosses_rule_sets", "rationale"])
            for prop in crep.proposals:
                w.writerow([prop.index, prop.name or "", prop.rule_set or "", prop.smarts, prop.action, ";".join(prop.evidence), ";".join(prop.covered_by),
                            ";".join(prop.covered_by_rule_sets), prop.crosses_rule_sets, prop.rationale])
        console.print("  Cleanup proposals: " + ", ".join(f"{v} {k}" for k, v in crep.by_action.items()) + f" (of {crep.n_patterns} patterns)")
        console.print(f"  [dim]{escape(crep.note)}[/dim]")
        console.print(f"[green]Saved:[/green] {cleanup}")
    if json_out:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(json.dumps(res.model_dump(), indent=2))
        console.print(f"[green]Full result saved to:[/green] {json_out}")


@app.command(name="smartsdiff")
def smartsdiff_cmd(
    old: str = typer.Argument(..., help="Catalogue A: pattern file, or rdkit:<NAME> (e.g. rdkit:PAINS)."),
    new: str = typer.Argument(..., help="Catalogue B: pattern file, or rdkit:<NAME>."),
    library: Optional[Path] = typer.Option(None, "--library", "-l", help="Reference molecules (.smi/.csv/.sdf). Default: bundled ChEMBL-derived set."),
    prep: str = typer.Option("implicit-h", "--prep", help="Molecule preparation: implicit-h, explicit-h or kekule."),
    max_molecules: Optional[int] = typer.Option(None, "--max-molecules", help="Use only the first N reference molecules."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Per-pattern table as CSV."),
    json_out: Optional[Path] = typer.Option(None, "--json", help="Full result as JSON."),
    show: int = typer.Option(15, "--show", help="Changed patterns to list in the terminal."),
) -> None:
    """Compare two SMARTS catalogues by what they do, not just what they say.

    Patterns are paired by name, then identical SMARTS, then identical hit set; each pair is
    classified as same hits / broadened / narrowed / shifted / broken / repaired, unpaired
    patterns as added / removed. The headline figure is the number of reference molecules whose
    flagged / not-flagged verdict differs between the two catalogues. Either side may be one of
    RDKit's built-in catalogues (rdkit:PAINS, rdkit:BRENK, rdkit:NIH, …), which are compared on
    hits alone because RDKit does not expose their SMARTS text.
    """
    import csv
    import json

    from chemlitmus.core.smartsaudit import PREPARATIONS
    from chemlitmus.core.smartsdiff import RDKIT_CATALOGS, diff_smarts

    if prep not in PREPARATIONS:
        console.print(f"[red]Error:[/red] Unknown --prep {prep!r}. Valid: {', '.join(PREPARATIONS)}.")
        raise typer.Exit(code=1)
    for side in (old, new):
        if side.lower().startswith("rdkit:"):
            if side.split(":", 1)[1] not in RDKIT_CATALOGS:
                console.print(f"[red]Error:[/red] Unknown RDKit catalogue {side!r}. Available: {', '.join('rdkit:' + c for c in RDKIT_CATALOGS)}.")
                raise typer.Exit(code=1)
        elif not Path(side).exists():
            console.print(f"[red]Error:[/red] Pattern file not found: {side}")
            raise typer.Exit(code=1)

    with console.status("[bold green]Matching both catalogues against the reference library...[/bold green]"):
        res = diff_smarts(old, new, library=library, preparation=prep, max_molecules=max_molecules)

    t1 = Table(title="[bold]SMARTS catalogue diff[/bold]", show_header=True, header_style="bold magenta")
    t1.add_column("", style="cyan"); t1.add_column("A", justify="right"); t1.add_column("B", justify="right")
    t1.add_row("Catalogue", escape(res.source_a), escape(res.source_b))
    t1.add_row("Patterns", str(res.n_patterns_a), str(res.n_patterns_b))
    t1.add_row("Molecules flagged", f"{res.flagged_a} ({res.flagged_a / max(res.n_molecules, 1):.1%})", f"{res.flagged_b} ({res.flagged_b / max(res.n_molecules, 1):.1%})")
    console.print(t1)
    console.print(f"  Reference: {res.n_molecules} molecules ({escape(res.library_source)}), preparation [bold]{res.preparation}[/bold]")
    console.print(f"  Paired {res.n_paired} patterns " + "(" + ", ".join(f"{v} by {k}" for k, v in res.paired_by.items()) + ")" if res.paired_by else f"  Paired {res.n_paired} patterns")
    if res.text_counts:
        console.print("  Text: " + ", ".join(f"{v} {k}" for k, v in res.text_counts.items()))
    console.print("  Behaviour: " + ", ".join(f"[bold]{v}[/bold] {k}" for k, v in sorted(res.semantic_counts.items(), key=lambda kv: -kv[1])))
    colour = "green" if res.verdict_changes == 0 else "yellow"
    console.print(f"  [{colour}]Verdict changes: {res.verdict_changes} molecules ({res.verdict_change_fraction:.1%})[/{colour}] — "
                  f"{res.flagged_only_b} newly flagged by B, {res.flagged_only_a} no longer flagged; {res.flagged_both} flagged by both")

    changed = [d for d in res.patterns if d.semantic_status != "same hits"]
    if changed:
        any_smarts = any((d.a and d.a.smarts) or (d.b and d.b.smarts) for d in changed[:show])
        t2 = Table(title=f"Changed patterns (first {min(show, len(changed))} of {len(changed)})", show_header=True, header_style="bold yellow")
        t2.add_column("Pattern", style="cyan", overflow="fold", min_width=18); t2.add_column("Change", min_width=9); t2.add_column("Text", min_width=9)
        t2.add_column("Hits A", justify="right"); t2.add_column("Hits B", justify="right"); t2.add_column("+", justify="right"); t2.add_column("−", justify="right")
        t2.add_column("Jaccard", justify="right")
        if any_smarts:
            t2.add_column("SMARTS A → B", overflow="fold", max_width=48)
        for d in changed[:show]:
            sm = ""
            if d.a and d.a.smarts and d.b and d.b.smarts and d.text_status == "rewritten":
                sm = f"{d.a.smarts} → {d.b.smarts}"
            elif (d.a and d.a.smarts) and not d.b:
                sm = d.a.smarts
            elif (d.b and d.b.smarts) and not d.a:
                sm = d.b.smarts
            row = [escape(d.key), d.semantic_status, d.text_status or "", str(d.hits_a), str(d.hits_b), str(d.gained), str(d.lost),
                   f"{d.jaccard:.2f}" if d.jaccard is not None else ""]
            t2.add_row(*(row + [escape(sm)] if any_smarts else row))
        console.print(t2)
    else:
        console.print("[green]Every paired pattern has the same hit set; nothing added or removed.[/green]")

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["key", "paired_by", "semantic_status", "text_status", "name_a", "smarts_a", "rule_set_a", "parses_a", "hits_a",
                        "name_b", "smarts_b", "rule_set_b", "parses_b", "hits_b", "gained", "lost", "jaccard", "examples_gained", "examples_lost"])
            for d in res.patterns:
                a, b = d.a, d.b
                w.writerow([d.key, d.paired_by or "", d.semantic_status, d.text_status or "",
                            a.name if a else "", a.smarts if a else "", a.rule_set if a else "", a.parses if a else "", d.hits_a,
                            b.name if b else "", b.smarts if b else "", b.rule_set if b else "", b.parses if b else "", d.hits_b,
                            d.gained, d.lost, "" if d.jaccard is None else f"{d.jaccard:.4f}", " ".join(d.examples_gained), " ".join(d.examples_lost)])
        console.print(f"[green]Saved:[/green] {output}")
    if json_out:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(json.dumps(res.model_dump(), indent=2))
        console.print(f"[green]Saved:[/green] {json_out}")
    raise typer.Exit(code=0)


@app.command(name="smartsproof")
def smartsproof_cmd(
    patterns: Optional[Path] = typer.Argument(None, help="Pattern file to prove redundancy within (CSV/TSV/XLSX with a 'smarts' column, or one SMARTS per line)."),
    subsumes_pair: Optional[Tuple[str, str]] = typer.Option(None, "--subsumes", help="A B: prove that every molecule matched by A is matched by B.", metavar="A B"),
    equivalent_pair: Optional[Tuple[str, str]] = typer.Option(None, "--equivalent", help="A B: prove that A and B match exactly the same molecules.", metavar="A B"),
    satisfiable_smarts: Optional[str] = typer.Option(None, "--satisfiable", help="Check whether any atom/bond state can satisfy every expression in the pattern."),
    max_container_atoms: int = typer.Option(40, "--max-container-atoms", help="Skip containers larger than this (reported as undecided)."),
    step_budget: int = typer.Option(200000, "--step-budget", help="Backtracking steps per pattern pair before giving up (undecided)."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Catalogue mode: per-pattern table as CSV."),
    json_out: Optional[Path] = typer.Option(None, "--json", help="Full result as JSON."),
    show: int = typer.Option(15, "--show", help="Catalogue mode: proven-redundant patterns to list."),
) -> None:
    """Prove SMARTS containment statically — no reference molecules needed.

    A proof is an injective atom mapping (the witness) under which every atom and bond expression
    of the broader pattern is implied by the one it maps onto; implication is decided exactly over
    a finite universe of atom states. A witness proves containment for every molecule. No witness
    proves nothing: the pair is undecided. Recursive SMARTS, isotopes, valence and ring-size
    primitives are outside the universe and reported as not analysable.
    """
    import csv
    import json

    from chemlitmus.core.smartsproof import equivalent, prove_catalogue, satisfiable, subsumes

    modes = sum(x is not None for x in (patterns, subsumes_pair, equivalent_pair, satisfiable_smarts))
    if modes != 1:
        console.print("[red]Error:[/red] Give exactly one of: a pattern file, --subsumes A B, --equivalent A B, --satisfiable S.")
        raise typer.Exit(code=1)

    def _witness_str(w, a_sm, b_sm):
        if not w:
            return ""
        return ", ".join(f"B[{k}]→A[{v}]" for k, v in sorted(w.items()))

    if subsumes_pair or equivalent_pair:
        a, b = subsumes_pair or equivalent_pair
        for sm in (a, b):
            if Chem.MolFromSmarts(sm) is None:
                console.print(f"[red]Error:[/red] Invalid SMARTS: {escape(sm)}")
                raise typer.Exit(code=1)
        r = subsumes(a, b) if subsumes_pair else equivalent(a, b)
        rel = "⊆" if subsumes_pair else "≡"
        colour = "green" if r.proven else ("yellow" if r.status.startswith("no witness") or r.status.startswith("budget") else "red")
        console.print(f"[bold]A[/bold] = {escape(a)}\n[bold]B[/bold] = {escape(b)}")
        console.print(f"[{colour}]A {rel} B: {r.status}[/{colour}]  [dim](atom universe: {r.universe_size:,} states)[/dim]")
        if r.witness:
            console.print(f"  Witness A ⊆ B: {_witness_str(r.witness, a, b)}  [dim](B atom → A atom)[/dim]")
        if equivalent_pair and r.witness_reverse:
            console.print(f"  Witness B ⊆ A: {', '.join(f'A[{k}]→B[{v}]' for k, v in sorted(r.witness_reverse.items()))}")
        if not r.proven and not r.status.startswith("not analysable"):
            console.print("  [dim]No witness is not a refutation: the containment may hold but is not provable by this method.[/dim]")
        if json_out:
            json_out.write_text(json.dumps(r.model_dump(), indent=2)); console.print(f"[green]Saved:[/green] {json_out}")
        raise typer.Exit(code=0 if r.proven else 1)

    if satisfiable_smarts:
        if Chem.MolFromSmarts(satisfiable_smarts) is None:
            console.print(f"[red]Error:[/red] Invalid SMARTS: {escape(satisfiable_smarts)}")
            raise typer.Exit(code=1)
        r = satisfiable(satisfiable_smarts)
        colour = {"satisfiable at atom and bond level": "green", "unsatisfiable": "red"}.get(r.status, "yellow")
        console.print(f"{escape(satisfiable_smarts)}: [{colour}]{r.status}[/{colour}]")
        if r.atom_index is not None:
            console.print(f"  Atom {r.atom_index} ({escape(r.atom_expression or '')}) matches no atom state.")
        if json_out:
            json_out.write_text(json.dumps(r.model_dump(), indent=2)); console.print(f"[green]Saved:[/green] {json_out}")
        raise typer.Exit(code=0 if r.satisfiable else 1)

    triples = load_patterns(patterns)
    if not triples:
        console.print("[red]Error:[/red] No patterns found.")
        raise typer.Exit(code=1)
    labels = [(sm, name or f"#{i}") for i, (sm, name, _) in enumerate(triples)]
    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"), BarColumn(), TaskProgressColumn(), console=console) as progress:
        task = progress.add_task("Proving...", total=len(labels))
        res = prove_catalogue(labels, max_container_atoms=max_container_atoms, step_budget=step_budget,
                              progress_callback=lambda n, total: progress.update(task, completed=n))

    n = res.n_patterns or 1
    t1 = Table(title="[bold]Static containment proofs[/bold]", show_header=True, header_style="bold magenta")
    t1.add_column("Outcome", style="cyan"); t1.add_column("Patterns", justify="right"); t1.add_column("Share", justify="right"); t1.add_column("Meaning")
    t1.add_row("Proven redundant", str(res.n_proven_redundant), f"{res.n_proven_redundant / n:.1%}", "Another catalogue pattern provably contains it (lower bound)")
    t1.add_row("  of which proven equivalent", str(res.n_proven_equivalent), f"{res.n_proven_equivalent / n:.1%}", "Containment proven both ways")
    t1.add_row("Undecided", str(res.n_no_witness), f"{res.n_no_witness / n:.1%}", "No witness found — not refuted")
    t1.add_row("Unsatisfiable", str(res.n_unsatisfiable), f"{res.n_unsatisfiable / n:.1%}", "An atom expression matches no atom state")
    t1.add_row("Not analysable", str(res.n_not_analysable), f"{res.n_not_analysable / n:.1%}", escape("; ".join(f"{k} ×{v}" for k, v in sorted(res.not_analysable_reasons.items(), key=lambda kv: -kv[1]))))
    console.print(t1)
    console.print(f"  {res.n_patterns} patterns · atom universe {res.universe_size:,} states · {res.n_pairs_tested:,} candidate pairs tested · {res.n_pairs_budget_exceeded} over budget")

    shown = [p for p in res.patterns if p.status in ("proven redundant", "proven equivalent")][:show]
    if shown:
        t2 = Table(title=f"Proven redundant (first {len(shown)})", show_header=True, header_style="bold yellow")
        t2.add_column("Pattern", style="cyan", overflow="fold", min_width=16); t2.add_column("SMARTS", overflow="fold", max_width=40)
        t2.add_column("Contained in", overflow="fold", min_width=16); t2.add_column("Witness", overflow="fold")
        for p in shown:
            t2.add_row(escape(p.name or ""), escape(p.smarts), escape(", ".join(p.proven_subsumed_by[:3]) + (" …" if len(p.proven_subsumed_by) > 3 else "")),
                       escape(", ".join(f"{k}→{v}" for k, v in sorted((p.witness or {}).items()))))
        console.print(t2)
        console.print("  [dim]Witness: container atom → pattern atom. Redundancy within one catalogue is a defect; across catalogues it is provenance.[/dim]")

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["index", "name", "smarts", "status", "reason", "proven_subsumed_by", "proven_equivalent_to", "witness", "n_undecided_pairs", "n_budget_exceeded"])
            for p in res.patterns:
                w.writerow([p.index, p.name or "", p.smarts, p.status, p.reason or "", ";".join(p.proven_subsumed_by), ";".join(p.proven_equivalent_to),
                            ";".join(f"{k}>{v}" for k, v in sorted((p.witness or {}).items())), p.n_undecided_pairs, p.n_budget_exceeded])
        console.print(f"[green]Saved:[/green] {output}")
    if json_out:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(json.dumps(res.model_dump(), indent=2)); console.print(f"[green]Saved:[/green] {json_out}")
    raise typer.Exit(code=0)
