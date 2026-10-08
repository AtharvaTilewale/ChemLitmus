"""Database lookups commands."""

from pathlib import Path
from typing import Optional
import typer
from rich.table import Table
from rich.markup import escape
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from chemlitmus import (
    lookup,
)
from chemlitmus.core.pubchem import PubChemCompound
from chemlitmus.utils.parsers import parse_compounds_file
from chemlitmus.utils.export import export_results
from chemlitmus.cli._app import app, console


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
    no_cache: bool = typer.Option(False, "--no-cache", help="Bypass the local record cache."),
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
        if res.disagreement:
            console.print(f"  [yellow]Differ by:[/yellow] {escape(res.disagreement)}  [dim](strictest shared identity level: {res.agreement_level})[/dim]")
            for p in res.pairwise:
                console.print(f"    [dim]{p.source_a} vs {p.source_b}:[/dim] {escape(p.description)}")
        for c in res.corrections:
            fix = f"replaced by {c.corrected_id}" if c.corrected else "no structure match; dropped"
            console.print(f"  [yellow]Corrected:[/yellow] {c.source} text search returned {escape(c.text_hit_id or '?')} "
                          f"({escape(c.text_hit_name or '')}) — {escape(c.difference)} — {fix}")
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
            res = resolve(query, sources=src_list, query_type=query_type, unichem=not no_unichem, timeout=timeout, use_cache=not no_cache)
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
            results.append(resolve(q, sources=src_list, query_type=query_type, unichem=not no_unichem, timeout=timeout, use_cache=not no_cache))
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

@app.command(name="concordance")
def concordance_cmd(
    file: Path = typer.Option(..., "--file", "-f", help="Queries, one per line (or CSV/TSV with a name/id/smiles column)."),
    sources: str = typer.Option("pubchem,chembl,chebi", "--sources", "-s", help="Comma-separated: pubchem, chembl, chebi, kegg, or all."),
    query_type: str = typer.Option("auto", "--type", "-t", help="auto, name, smiles, inchikey, or id."),
    timeout: float = typer.Option(40.0, "--timeout", help="Seconds allowed per query."),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Per-query table as CSV."),
    json_out: Optional[Path] = typer.Option(None, "--json", help="Full report as JSON."),
    show: int = typer.Option(15, "--show", help="Disagreeing queries to list in the terminal."),
    no_cache: bool = typer.Option(False, "--no-cache", help="Bypass the local record cache."),
) -> None:
    """Measure how often databases return the same structure for the same name.

    Resolves every query across the chosen sources and reports, per query, whether the sources
    agree (by InChIKey), the strictest identity level they share, and what differs — salt or
    charge form, tautomer, stereochemistry, or a different compound. The summary also counts how
    often each database's text search landed on a different structure than the majority and
    had to be corrected through the consensus InChIKey.
    """
    import csv
    import json

    from chemlitmus.providers import PROVIDERS, concordance

    src_list = [s.strip().lower() for s in sources.split(",") if s.strip()]
    bad = [s for s in src_list if s != "all" and s not in PROVIDERS]
    if bad:
        console.print(f"[red]Error:[/red] Unknown source(s) {bad}. Available: {sorted(PROVIDERS)} or 'all'.")
        raise typer.Exit(code=1)
    if query_type.lower() not in ("auto", "name", "smiles", "inchikey", "id"):
        console.print(f"[red]Error:[/red] Unknown --type {query_type!r}. Valid: auto, name, smiles, inchikey, id.")
        raise typer.Exit(code=1)
    queries = parse_compounds_file(file)
    if not queries:
        console.print("[red]Error:[/red] No queries found in the input file.")
        raise typer.Exit(code=1)

    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"), BarColumn(), TaskProgressColumn(), console=console) as progress:
        task = progress.add_task("Resolving...", total=len(queries))
        rep = concordance(queries, sources=src_list, query_type=query_type, timeout=timeout, use_cache=not no_cache,
                          progress_callback=lambda n, total: progress.update(task, completed=n))

    t1 = Table(title=f"[bold]Concordance across {', '.join(rep.sources)}[/bold]", show_header=True, header_style="bold magenta")
    t1.add_column("Measure", style="cyan"); t1.add_column("Count", justify="right"); t1.add_column("Share", justify="right")
    n = rep.n_queries or 1
    t1.add_row("Queries", str(rep.n_queries), "")
    t1.add_row("Found in at least one source", str(rep.n_found_any), f"{rep.n_found_any / n:.1%}")
    t1.add_row("Found in every source", str(rep.n_found_all), f"{rep.n_found_all / n:.1%}")
    for k in ("agree", "disagree", "unknown", "not found"):
        if rep.agreement_counts.get(k):
            t1.add_row(f"Agreement: {k}", str(rep.agreement_counts[k]), f"{rep.agreement_counts[k] / n:.1%}")
    console.print(t1)

    if rep.level_counts:
        t2 = Table(title="Strictest identity level shared by all returned structures", show_header=True, header_style="bold blue")
        t2.add_column("Level", style="cyan"); t2.add_column("Queries", justify="right"); t2.add_column("Meaning")
        meaning = {"exact": "identical structures", "parent": "differ by salt, counter-ion or charge form", "tautomer": "differ by tautomer",
                   "nostereo": "differ by stereochemistry", "skeleton": "differ by stereochemistry and tautomer",
                   "formula": "same formula, different constitution", "different compounds": "no shared level"}
        for lv in ("exact", "parent", "tautomer", "nostereo", "skeleton", "formula", "different compounds"):
            if rep.level_counts.get(lv):
                t2.add_row(lv, str(rep.level_counts[lv]), meaning[lv])
        console.print(t2)

    t3 = Table(title="Per source", show_header=True, header_style="bold blue")
    t3.add_column("Source", style="cyan"); t3.add_column("Found", justify="right"); t3.add_column("Errors", justify="right")
    t3.add_column("Text hit corrected by structure", justify="right"); t3.add_column("What the text hit was", overflow="fold")
    for s in rep.sources:
        cls = rep.correction_classes_by_source.get(s, {})
        t3.add_row(s, str(rep.found_by_source.get(s, 0)), str(rep.errors_by_source.get(s, 0)), str(rep.corrections_by_source.get(s, 0)),
                   escape("; ".join(f"{k} ×{v}" for k, v in sorted(cls.items(), key=lambda kv: -kv[1]))))
    console.print(t3)

    shown = [r for r in rep.rows if r.agreement == "disagree" or r.corrections][:show]
    if shown:
        t4 = Table(title=f"Disagreements and corrections (first {len(shown)})", show_header=True, header_style="bold yellow",
                   caption="a → b: text search returned a, replaced by b (the record for the consensus structure); a → ∅: no such record, a kept")
        t4.add_column("Query", style="cyan", min_width=12, overflow="fold"); t4.add_column("Level", min_width=8, overflow="fold")
        t4.add_column("Differs by", min_width=14, overflow="fold")
        for s in rep.sources:
            t4.add_column(s, min_width=14, overflow="fold")
        for r in shown:
            cells = []
            for s in rep.sources:
                cid = r.ids.get(s, r.errors.get(s, ""))
                corr = next((c for c in r.corrections if c.source == s), None)
                if corr:
                    cid = f"{corr.text_hit_id} → {corr.corrected_id or '∅'}"
                cells.append(escape(str(cid)))
            t4.add_row(escape(r.query), r.agreement_level or "—", escape(r.disagreement or ""), *cells)
        console.print(t4)

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["query", "n_found", "agreement", "agreement_level", "disagreement", "consensus_inchikey"]
                       + [f"{s}_id" for s in rep.sources] + [f"{s}_name" for s in rep.sources] + [f"{s}_inchikey" for s in rep.sources]
                       + [f"{s}_text_hit" for s in rep.sources] + [f"{s}_correction" for s in rep.sources] + [f"{s}_status" for s in rep.sources])
            for r in rep.rows:
                corr = {c.source: c for c in r.corrections}
                w.writerow([r.query, r.n_found, r.agreement, r.agreement_level or "", r.disagreement or "", r.consensus_inchikey or ""]
                           + [r.ids.get(s, "") for s in rep.sources] + [r.names.get(s, "") for s in rep.sources] + [r.inchikeys.get(s, "") for s in rep.sources]
                           + [corr[s].text_hit_id if s in corr else "" for s in rep.sources]
                           + [corr[s].difference if s in corr else "" for s in rep.sources]
                           + [("found" if s in r.ids else r.errors.get(s, "")) for s in rep.sources])
        console.print(f"[green]Saved:[/green] {output}")
    if json_out:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(json.dumps(rep.model_dump(), indent=2, default=str))
        console.print(f"[green]Saved:[/green] {json_out}")


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
