"""Screening, similarity and substructure search commands."""

from pathlib import Path
from typing import Optional
import typer
from rich.table import Table
from rich.markup import escape
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from chemlitmus import (
    compute_fingerprint, apply_filters, compute_similarity,
    substructure_search, extract_scaffold, rgroup_decomposition,
)
from chemlitmus.logging_config import logger
from chemlitmus.cli._app import app, console


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
