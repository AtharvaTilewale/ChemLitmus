"""Dataset audit, leakage and label-conflict commands."""

from pathlib import Path
from typing import List, Optional

import typer
from rich.markup import escape
from rich.progress import BarColumn, Progress, SpinnerColumn, TaskProgressColumn, TextColumn
from rich.table import Table

from chemlitmus.cli._app import app, console
from chemlitmus.core.records import SchemaError


def _roles_from(structure_column, id_column, split_column, endpoint_column, units_column, relation_column, date_column, source_column, target_column):
    roles = {}
    for role, col in (("split", split_column), ("endpoint", endpoint_column), ("units", units_column), ("relation", relation_column),
                      ("date", date_column), ("source", source_column), ("target", target_column)):
        if col:
            roles[role] = col
    return roles


@app.command(name="audit")
def audit_cmd(
    dataset: Path = typer.Argument(..., help="CSV/TSV/XLSX/SMI/SDF dataset."),
    output_dir: Path = typer.Option(Path("audit_results"), "--output-dir", "-o", help="Where to write records.csv, issues.csv, identity_groups.csv, summary.json, audit.json, policy.json, manifest.json, report.html."),
    config: Optional[str] = typer.Option(None, "--config", "-c", help="Chemical policy: 'parent' (default), 'conservative', or a policy JSON file."),
    structure_column: Optional[str] = typer.Option(None, "--structure-column", help="Column holding SMILES (default: detected from header)."),
    id_column: Optional[str] = typer.Option(None, "--id-column", help="Column holding the source id."),
    split_column: Optional[str] = typer.Option(None, "--split-column", help="Column with split labels; enables leakage checks."),
    endpoint_column: Optional[str] = typer.Option(None, "--endpoint-column", help="Label / measurement column; enables conflict checks."),
    units_column: Optional[str] = typer.Option(None, "--units-column"),
    relation_column: Optional[str] = typer.Option(None, "--relation-column", help="Column with =, <, >, <=, >= qualifiers."),
    date_column: Optional[str] = typer.Option(None, "--date-column", help="Enables the temporal-order check between splits."),
    source_column: Optional[str] = typer.Option(None, "--source-column", help="Endpoint context: measurements are compared only within the same source."),
    target_column: Optional[str] = typer.Option(None, "--target-column", help="Endpoint context: compared only within the same target."),
    first_column: bool = typer.Option(False, "--first-column", help="Use the first column as the structure when no header is recognised (otherwise a schema error)."),
    clean: bool = typer.Option(False, "--clean", help="Also write clean.csv (standardised structures, exclusions listed in exclusions.csv)."),
    exclude_codes: str = typer.Option("DUP_EXACT", "--exclude-codes", help="Clean export: comma-separated issue codes that exclude a record."),
    exclude_severity: Optional[str] = typer.Option(None, "--exclude-severity", help="Clean export: exclude records whose max severity is at least this (warning|error)."),
    fail_on_severity: Optional[str] = typer.Option(None, "--fail-on-severity", help="Gate: exit 3 when any issue is at least this severity (warning|error)."),
    max_invalid_fraction: Optional[float] = typer.Option(None, "--max-invalid-fraction", help="Gate: exit 3 when (empty+invalid+error)/total exceeds this."),
    max_duplicate_fraction: Optional[float] = typer.Option(None, "--max-duplicate-fraction", help="Gate: exit 3 when collapsed/annotated exceeds this."),
    no_split_overlap: bool = typer.Option(False, "--no-split-overlap", help="Gate: exit 3 on any split overlap at the policy identity level."),
    no_label_conflicts: bool = typer.Option(False, "--no-label-conflicts", help="Gate: exit 3 on any label conflict."),
    require_complete: bool = typer.Option(False, "--require-complete", help="Gate: exit 4 when any computation failed."),
) -> None:
    """Audit a chemical dataset: parsing, standardisation, identity groups, structural flags,
    descriptors, alerts, split leakage and label conflicts — with every record accounted for.

    Exit codes: 0 pass · 1 configuration error · 3 gate violated · 4 processing incomplete (with --require-complete).
    """
    from chemlitmus.core.dataset_audit import audit_dataset
    from chemlitmus.core.policy import resolve_policy
    from chemlitmus.core.report import CleanPolicy, GatePolicy, evaluate_gates, write_audit_outputs

    for sev, name in ((fail_on_severity, "--fail-on-severity"), (exclude_severity, "--exclude-severity")):
        if sev and sev not in ("warning", "error"):
            console.print(f"[red]Error:[/red] {name} must be warning or error.")
            raise typer.Exit(code=1)
    if not dataset.exists():
        console.print(f"[red]Error:[/red] File not found: {dataset}")
        raise typer.Exit(code=1)
    try:
        policy = resolve_policy(config)
    except Exception as exc:
        console.print(f"[red]Error:[/red] Invalid policy: {escape(str(exc))}")
        raise typer.Exit(code=1)
    roles = _roles_from(structure_column, id_column, split_column, endpoint_column, units_column, relation_column, date_column, source_column, target_column)
    try:
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"), BarColumn(), TaskProgressColumn(), console=console) as progress:
            task = progress.add_task("Auditing...", total=None)
            audit = audit_dataset(str(dataset), policy=policy, structure_column=structure_column, id_column=id_column, roles=roles or None,
                                  ambiguous="first" if first_column else "error",
                                  progress_callback=lambda n, total: progress.update(task, total=total, completed=n))
    except (SchemaError, ValueError) as exc:
        console.print(f"[red]Error:[/red] {escape(str(exc))}")
        raise typer.Exit(code=1)
    gates = GatePolicy(fail_on_severity=fail_on_severity, max_invalid_fraction=max_invalid_fraction, max_duplicate_fraction=max_duplicate_fraction,
                       no_split_overlap=no_split_overlap, no_label_conflicts=no_label_conflicts, require_complete_processing=require_complete)
    clean_pol = CleanPolicy(exclude_codes=[c.strip() for c in exclude_codes.split(",") if c.strip()], exclude_min_severity=exclude_severity) if clean else None
    paths = write_audit_outputs(audit, output_dir, input_path=dataset, gates=gates, clean=clean_pol, command="chemlitmus audit " + " ".join(escape(a) for a in [str(dataset)]))
    gate = evaluate_gates(audit, gates)
    _print_audit(audit, gate)
    console.print(f"[green]Outputs:[/green] {output_dir}/  (" + ", ".join(Path(p).name for p in paths.values()) + ")")
    raise typer.Exit(code=gate.exit_code)


def _print_audit(audit, gate) -> None:
    s = audit.summary
    t = Table(title=f"[bold]Dataset audit — {escape(Path(audit.source_file).name)}[/bold]", show_header=True, header_style="bold magenta")
    t.add_column("Measure", style="cyan"); t.add_column("Value", justify="right"); t.add_column("Note")
    t.add_row("Records", str(s.n_total), f"ok {s.n_ok} · empty {s.n_empty} · invalid {s.n_invalid} · unsupported {s.n_unsupported} · error {s.n_error}")
    t.add_row("Annotated", str(s.n_annotated), "parsed, standardised, identity computed")
    t.add_row("Processing complete", str(s.processing_complete), "" if s.processing_complete else "some computations failed — not negatives")
    t.add_row(f"Distinct compounds ({s.identity_level})", str(s.n_groups), f"{s.n_collapsed} records collapse onto an earlier one")
    for sev in ("error", "warning", "info"):
        if s.issues_by_severity.get(sev):
            t.add_row(f"Issues: {sev}", str(s.issues_by_severity[sev]), ", ".join(f"{c} ×{n}" for c, n in sorted(s.issues_by_code.items(), key=lambda kv: -kv[1]) if _sev_of(c) == sev)[:90])
    if s.alerts_by_set:
        t.add_row("Alert matches", ", ".join(f"{k} {v}" for k, v in s.alerts_by_set.items()), f"of {s.alert_denominator} annotated; context, not a verdict")
    if audit.leakage:
        L = audit.leakage
        for p in L["pairs"]:
            ov = {o["level"]: o["n_eval_records"] for o in p["overlap"]}
            t.add_row(f"Leakage {p['reference_split']}→{p['evaluation_split']}", str(ov.get(L["identity_level"], 0)),
                      f"exact {ov['exact']} · parent {ov['parent']} · tautomer {ov['tautomer']} · nostereo {ov['nostereo']} · skeleton {ov['skeleton']} (nested, not additive); related {p['n_related_by_similarity']} by similarity, {p['scaffold_overlap_records']} by scaffold")
    if audit.label_conflicts:
        C = audit.label_conflicts
        t.add_row("Label conflicts", str(C["n_conflicts"]), f"{C['n_groups_compared']} groups compared ({C['kind']}); censored {C['n_censored']}, missing units {C['n_missing_units']}")
    console.print(t)
    for w in s.dataset_warnings:
        console.print(f"  [dim]{escape(w)}[/dim]")
    colour = "green" if gate.passed else "red"
    console.print(f"[{colour}]Gate: {'PASS' if gate.passed else 'FAIL'}[/{colour}]" + (" — " + "; ".join(escape(v) for v in gate.violations) if gate.violations else ""))


def _sev_of(code: str) -> str:
    from chemlitmus.core.dataset_audit import ISSUE_CATALOGUE
    return ISSUE_CATALOGUE.get(code, ("info",))[0]


@app.command(name="leakage")
def leakage_cmd(
    files: List[Path] = typer.Argument(None, help="Two or more split files (first = reference/training), or one file with --split-column."),
    split_column: Optional[str] = typer.Option(None, "--split-column", help="Single-file mode: column with split labels."),
    reference: Optional[str] = typer.Option(None, "--reference", help="Split treated as training (default: the largest)."),
    structure_column: Optional[str] = typer.Option(None, "--structure-column"),
    id_column: Optional[str] = typer.Option(None, "--id-column"),
    date_column: Optional[str] = typer.Option(None, "--date-column"),
    config: Optional[str] = typer.Option(None, "--config", "-c", help="Chemical policy (identity level, fingerprint, similarity threshold)."),
    json_out: Optional[Path] = typer.Option(None, "--json"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="CSV of overlapping / related evaluation records."),
    fail_on_overlap: bool = typer.Option(False, "--fail-on-overlap", help="Exit 3 on any overlap at the policy identity level."),
) -> None:
    """Report chemical leakage between splits as distinct, non-additive evidence classes
    (exact, parent, tautomer, nostereo, skeleton), plus scaffold and similarity relatedness."""
    import csv
    import json

    from chemlitmus.core.leakage import leakage_report
    from chemlitmus.core.policy import resolve_policy
    from chemlitmus.core.records import read_records

    policy = resolve_policy(config)
    splits = {}
    try:
        if split_column:
            if not files or len(files) != 1:
                console.print("[red]Error:[/red] --split-column needs exactly one file."); raise typer.Exit(code=1)
            rs = read_records(files[0], structure_column=structure_column, id_column=id_column, roles={"split": split_column, **({"date": date_column} if date_column else {})})
            for r in rs.ok_records():
                lab = str(r.fields.get(split_column, "")).strip()
                if lab:
                    splits.setdefault(lab, []).append((r.source_id or r.record_id, r.parsed_smiles, str(r.fields.get(date_column)).strip() if date_column and r.fields.get(date_column) else None))
            n_unassigned = sum(1 for r in rs.ok_records() if not str(r.fields.get(split_column, "")).strip())
            n_excluded = rs.n_total - rs.n_ok
        else:
            if not files or len(files) < 2:
                console.print("[red]Error:[/red] Give two or more files, or one file with --split-column."); raise typer.Exit(code=1)
            n_unassigned = n_excluded = 0
            for f in files:
                rs = read_records(f, structure_column=structure_column, id_column=id_column, roles={"date": date_column} if date_column else None)
                splits[f.stem] = [(r.source_id or r.record_id, r.parsed_smiles, str(r.fields.get(date_column)).strip() if date_column and r.fields.get(date_column) else None) for r in rs.ok_records()]
                n_excluded += rs.n_total - rs.n_ok
            reference = reference or files[0].stem
    except (SchemaError, FileNotFoundError, ValueError) as exc:
        console.print(f"[red]Error:[/red] {escape(str(exc))}"); raise typer.Exit(code=1)
    rep = leakage_report(splits, policy, reference=reference)
    rep.n_unassigned = n_unassigned
    console.print(f"[bold]Leakage report[/bold] · identity level [bold]{rep.identity_level}[/bold] · " + ", ".join(f"{k} {v}" for k, v in rep.n_by_split.items())
                  + f" · unassigned {rep.n_unassigned} · excluded (unparseable/empty) {n_excluded}")
    console.print("  within-split duplicates: " + ", ".join(f"{k} {v}" for k, v in rep.within_split_duplicates.items()))
    for p in rep.pairs:
        t = Table(title=f"{escape(p.reference_split)} → {escape(p.evaluation_split)}  ({p.n_reference} vs {p.n_evaluation})", show_header=True, header_style="bold blue")
        t.add_column("Evidence class"); t.add_column("Eval records", justify="right"); t.add_column("Shared keys", justify="right"); t.add_column("Fraction", justify="right")
        for o in p.overlap:
            mark = " [bold](policy level)[/bold]" if o.level == rep.identity_level else ""
            t.add_row(o.level + mark, str(o.n_eval_records), str(o.n_eval_groups), f"{o.fraction_of_eval:.1%}")
        t.add_row("[dim]formula match (not identity)[/dim]", str(p.formula_matches), "", "")
        t.add_row("[dim]scaffold overlap (relatedness)[/dim]", str(p.scaffold_overlap_records), f"{p.n_evaluation - p.n_acyclic_evaluation} cyclic", f"{p.scaffold_overlap_fraction:.1%}")
        t.add_row(f"[dim]similarity ≥ {p.neighbour_threshold} ({p.neighbour_fingerprint})[/dim]", str(p.n_related_by_similarity), "", escape(str(p.nearest_neighbour_similarity)))
        console.print(t)
        console.print(f"  [dim]temporal: {p.temporal_violations} violation(s) — {escape(p.temporal_note or '')}[/dim]")
    console.print(f"  [dim]{escape(rep.note)}[/dim]")
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh); w.writerow(["reference_split", "evaluation_split", "eval_record", "class", "reference_record", "detail"])
            for p in rep.pairs:
                for o in p.overlap:
                    for ex in o.examples:
                        w.writerow([p.reference_split, p.evaluation_split, ex["eval_record_id"], o.level, ex["reference_record_id"], ex["key"]])
                for n in p.neighbours:
                    if n.similarity >= p.neighbour_threshold:
                        w.writerow([p.reference_split, p.evaluation_split, n.eval_record_id, "similarity", n.reference_record_id, n.similarity])
        console.print(f"[green]Saved:[/green] {output}")
    if json_out:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(json.dumps(rep.model_dump(), indent=2)); console.print(f"[green]Saved:[/green] {json_out}")
    n_over = sum(o.n_eval_records for p in rep.pairs for o in p.overlap if o.level == rep.identity_level)
    if fail_on_overlap and n_over:
        console.print(f"[red]Gate: FAIL — {n_over} evaluation record(s) overlap at level {rep.identity_level}[/red]")
        raise typer.Exit(code=3)
    raise typer.Exit(code=0)


@app.command(name="conflicts")
def conflicts_cmd(
    dataset: Path = typer.Argument(..., help="Dataset with an endpoint column."),
    endpoint_column: str = typer.Option(..., "--endpoint-column", "-e"),
    units_column: Optional[str] = typer.Option(None, "--units-column"),
    relation_column: Optional[str] = typer.Option(None, "--relation-column"),
    context_columns: Optional[str] = typer.Option(None, "--context", help="Comma-separated columns (target, assay, source) that must match for measurements to be compared."),
    structure_column: Optional[str] = typer.Option(None, "--structure-column"),
    id_column: Optional[str] = typer.Option(None, "--id-column"),
    tolerance: float = typer.Option(1.0, "--tolerance", help="log10 units for molar quantities (1.0 = ten-fold), absolute otherwise."),
    kind: str = typer.Option("auto", "--kind", help="auto | classification | quantitative"),
    config: Optional[str] = typer.Option(None, "--config", "-c"),
    json_out: Optional[Path] = typer.Option(None, "--json"),
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="CSV of conflicting groups' measurements."),
) -> None:
    """Find contradictory labels or measurements among records of the same compound."""
    import csv
    import json

    from chemlitmus.core.identity import compute_identity
    from chemlitmus.core.labels import label_conflicts
    from chemlitmus.core.policy import resolve_policy
    from chemlitmus.core.records import read_records

    if kind not in ("auto", "classification", "quantitative"):
        console.print("[red]Error:[/red] --kind must be auto, classification or quantitative."); raise typer.Exit(code=1)
    policy = resolve_policy(config)
    ctx = [c.strip() for c in (context_columns or "").split(",") if c.strip()]
    try:
        roles = {"endpoint": endpoint_column, **({"units": units_column} if units_column else {}), **({"relation": relation_column} if relation_column else {})}
        rs = read_records(dataset, structure_column=structure_column, id_column=id_column, roles=roles)
        for c in ctx:
            if c not in rs.columns:
                raise SchemaError(f"Context column {c!r} not found. Columns: {rs.columns}")
    except (SchemaError, FileNotFoundError, ValueError) as exc:
        console.print(f"[red]Error:[/red] {escape(str(exc))}"); raise typer.Exit(code=1)
    recs = []
    for r in rs.ok_records():
        k = compute_identity(r.parsed_smiles)
        d = dict(r.fields); d.update({"record_id": r.record_id, "source_id": r.source_id, "key": k.key(policy.identity_level) if k.is_valid else None})
        recs.append(d)
    rep = label_conflicts(recs, policy, endpoint_field=endpoint_column, units_field=units_column, relation_field=relation_column, context_fields=ctx, tolerance=tolerance, kind=kind)
    console.print(f"[bold]Label conflicts[/bold] · endpoint [cyan]{escape(endpoint_column)}[/cyan] ({rep.kind}) · identity level {rep.identity_level} · context {escape(', '.join(ctx) or 'none')} · tolerance {tolerance}")
    console.print(f"  {rep.n_groups_compared} groups compared · [bold]{rep.n_conflicts} in conflict[/bold] ({rep.n_records_in_conflict} records) · censored {rep.n_censored} · missing units {rep.n_missing_units} · unparseable {rep.n_unparseable} · {rs.n_total - rs.n_ok} records not parsed")
    if rep.groups:
        t = Table(show_header=True, header_style="bold yellow")
        t.add_column("Group"); t.add_column("Records", justify="right"); t.add_column("Spread / labels"); t.add_column("Scale"); t.add_column("Measurements", overflow="fold")
        for g in rep.groups[:30]:
            ms = "; ".join(f"{m.source_id or m.record_id}: {m.relation if m.relation != '=' else ''}{m.raw_value} {m.units or ''}".strip() for m in g.measurements) if g.measurements else str(g.labels)
            t.add_row(g.group_id, str(g.n_records), str(g.spread if g.spread is not None else g.labels), g.scale or "", escape(ms))
        console.print(t)
    console.print(f"  [dim]{escape(rep.note)}[/dim]")
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh); w.writerow(["group_id", "key", "context", "record_id", "source_id", "relation", "raw_value", "units", "value_nm", "log_value", "censored", "comparable", "conversion", "labels"])
            for g in rep.groups:
                if g.measurements:
                    for m in g.measurements:
                        w.writerow([g.group_id, g.key, json.dumps(g.context), m.record_id, m.source_id or "", m.relation, m.raw_value, m.units or "", m.value_nm if m.value_nm is not None else "", m.log_value if m.log_value is not None else "", m.censored, m.comparable, m.conversion or "", ""])
                else:
                    w.writerow([g.group_id, g.key, json.dumps(g.context), "", "", "", "", "", "", "", "", "", "", json.dumps(g.labels)])
        console.print(f"[green]Saved:[/green] {output}")
    if json_out:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(json.dumps(rep.model_dump(), indent=2)); console.print(f"[green]Saved:[/green] {json_out}")
    raise typer.Exit(code=0)
