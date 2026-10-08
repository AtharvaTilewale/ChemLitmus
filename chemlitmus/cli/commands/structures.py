"""Structure retrieval and generation commands."""

from pathlib import Path
from typing import Optional
import re
import hashlib
import typer
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn

from chemlitmus import (
    lookup, download_structure, generate_structure, validate_smiles,
    validate_reaction, generate_conformers, augment_smiles, map_atoms,
)
from chemlitmus.utils.parsers import parse_compounds_file
from chemlitmus.cli._app import app, console


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
