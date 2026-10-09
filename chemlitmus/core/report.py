"""Audit outputs: per-record CSV, issues CSV, identity groups, summary JSON, authoritative JSON,
manifest, optional clean export with exclusions, and an offline HTML report."""

from __future__ import annotations

import base64
import csv
import html
import io
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from pydantic import BaseModel, Field

from chemlitmus.core.dataset_audit import ISSUE_CATALOGUE, SEVERITY_RANK, DatasetAudit
from chemlitmus.core.manifest import new_manifest

try:
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
except ImportError:  # pragma: no cover
    Chem = None

OUTPUT_SCHEMA_VERSION = "1"


# ------------------------------------------------------------------------------------ gates

class GatePolicy(BaseModel):
    fail_on_severity: Optional[str] = Field(None, description="Fail when any issue has at least this severity (warning | error).")
    max_invalid_fraction: Optional[float] = Field(None, description="Fail when (empty+invalid+error)/total exceeds this.")
    max_duplicate_fraction: Optional[float] = Field(None, description="Fail when collapsed records / annotated exceeds this (policy identity level).")
    no_split_overlap: bool = Field(False, description="Fail on any SPLIT_OVERLAP issue.")
    no_label_conflicts: bool = Field(False, description="Fail on any LABEL_CONFLICT issue.")
    require_complete_processing: bool = Field(False, description="Fail (exit 4) when processing_complete is false.")


class GateResult(BaseModel):
    passed: bool
    exit_code: int = Field(description="0 pass · 3 policy violation · 4 partial processing")
    violations: List[str] = Field(default_factory=list)
    partial_processing: bool = False


def evaluate_gates(audit: DatasetAudit, gates: GatePolicy) -> GateResult:
    s = audit.summary
    v: List[str] = []
    if gates.fail_on_severity:
        thr = SEVERITY_RANK[gates.fail_on_severity]
        n = sum(c for sev, c in s.issues_by_severity.items() if SEVERITY_RANK[sev] >= thr)
        if n:
            v.append(f"{n} issue(s) at severity >= {gates.fail_on_severity}")
    if gates.max_invalid_fraction is not None and s.n_total:
        frac = (s.n_empty + s.n_invalid + s.n_error) / s.n_total
        if frac > gates.max_invalid_fraction:
            v.append(f"invalid fraction {frac:.3f} > {gates.max_invalid_fraction}")
    if gates.max_duplicate_fraction is not None and s.n_annotated:
        frac = s.n_collapsed / s.n_annotated
        if frac > gates.max_duplicate_fraction:
            v.append(f"duplicate fraction {frac:.3f} at level {s.identity_level} > {gates.max_duplicate_fraction}")
    if gates.no_split_overlap and s.issues_by_code.get("SPLIT_OVERLAP"):
        v.append(f"{s.issues_by_code['SPLIT_OVERLAP']} split overlap(s) at level {s.identity_level}")
    if gates.no_label_conflicts and s.issues_by_code.get("LABEL_CONFLICT"):
        v.append(f"{s.issues_by_code['LABEL_CONFLICT']} record(s) in label conflict")
    partial = not s.processing_complete
    if v:
        return GateResult(passed=False, exit_code=3, violations=v, partial_processing=partial)
    if partial and gates.require_complete_processing:
        return GateResult(passed=False, exit_code=4, violations=["processing incomplete (see STD_FAILED / unevaluated alerts)"], partial_processing=True)
    return GateResult(passed=True, exit_code=0, partial_processing=partial)


# ------------------------------------------------------------------------------------ writers

def _flat(v) -> str:
    if v is None:
        return ""
    if isinstance(v, (list, dict)):
        return json.dumps(v, default=str)
    return str(v)


def write_records_csv(audit: DatasetAudit, path: Path) -> None:
    orig_cols: List[str] = []
    for a in audit.records:
        for k in a.fields:
            if k not in orig_cols:
                orig_cols.append(k)
    ann_cols = ["record_id", "position", "source_id", "status", "max_severity", "input_structure", "parsed_smiles", "standardized_smiles",
                "standardization_changed_identity", "repair_candidate", "repair_status", "group_id", "identity_exact", "identity_parent",
                "identity_tautomer", "identity_nostereo", "identity_skeleton", "formula", "n_components", "n_organic_components", "net_charge",
                "has_isotopes", "n_stereocentres", "n_unassigned_stereocentres", "scaffold", "split", "alerts", "issue_codes", "issues"]
    desc_cols = ["mw", "logp", "hbd", "hba", "tpsa", "rotatable_bonds", "heavy_atoms", "rings", "fraction_csp3"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([f"in:{c}" for c in orig_cols] + ann_cols + desc_cols)
        for a in audit.records:
            k = a.identity
            row = [_flat(a.fields.get(c)) for c in orig_cols] + [
                a.record_id, a.position, _flat(a.source_id), a.status, a.max_severity, _flat(a.input_structure), _flat(a.parsed_smiles), _flat(a.standardized_smiles),
                a.standardization_changed_identity, _flat(a.repair_candidate), _flat(a.repair_status), _flat(a.group_id),
                _flat(k.exact if k else None), _flat(k.parent if k else None), _flat(k.tautomer if k else None), _flat(k.nostereo if k else None), _flat(k.skeleton if k else None), _flat(k.formula if k else None),
                a.n_components, a.n_organic_components, _flat(a.net_charge), a.has_isotopes, a.n_stereocentres, a.n_unassigned_stereocentres, _flat(a.scaffold), _flat(a.split),
                ";".join(f"{x['set']}:{x['description']}" for x in a.alerts), ";".join(i.code for i in a.issues), " | ".join(f"[{i.severity}] {i.code}: {i.message}" for i in a.issues),
            ] + [_flat(a.descriptors.get(c)) for c in desc_cols]
            w.writerow(row)


def write_issues_csv(audit: DatasetAudit, path: Path) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["code", "severity", "record_id", "source_id", "message", "suggested_action", "evidence"])
        for i in sorted(audit.issues, key=lambda x: (-SEVERITY_RANK[x.severity], x.code, x.record_id or "")):
            w.writerow([i.code, i.severity, _flat(i.record_id), _flat(i.source_id), i.message, i.suggested_action, json.dumps(i.evidence, default=str)])


def write_groups_csv(audit: DatasetAudit, path: Path) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["group_id", "level", "size", "key", "record_ids", "source_ids", "differs_by", "splits"])
        for g in audit.groups:
            w.writerow([g.group_id, g.level, g.size, g.key, ";".join(g.record_ids), ";".join(_flat(s) for s in g.source_ids), ";".join(g.differs_by), ";".join(g.splits)])


class CleanPolicy(BaseModel):
    """What a 'clean' export includes. Every exclusion is listed with its reason."""

    exclude_statuses: List[str] = Field(default_factory=lambda: ["empty", "invalid", "unsupported", "error"])
    exclude_codes: List[str] = Field(default_factory=lambda: ["DUP_EXACT"], description="Records carrying any of these issue codes are excluded.")
    exclude_min_severity: Optional[str] = Field(None, description="Exclude records whose max severity is at least this.")
    use_standardized_structure: bool = Field(True, description="Write the standardised SMILES as the structure; the original is kept in its own column.")
    apply_repairs: bool = Field(False, description="Write candidate repairs for invalid records (requires repair.mode=apply in the chemical policy).")


def write_clean(audit: DatasetAudit, out_path: Path, excl_path: Path, clean: CleanPolicy) -> Tuple[int, int]:
    """Write the clean dataset and the exclusion table. Returns (n_written, n_excluded)."""
    orig_cols: List[str] = []
    for a in audit.records:
        for k in a.fields:
            if k not in orig_cols:
                orig_cols.append(k)
    n_w = n_x = 0
    with open(out_path, "w", newline="", encoding="utf-8") as fo, open(excl_path, "w", newline="", encoding="utf-8") as fx:
        wo, wx = csv.writer(fo), csv.writer(fx)
        wo.writerow(["record_id", "source_id", "structure", "original_structure", "structure_source", "group_id"] + orig_cols)
        wx.writerow(["record_id", "source_id", "position", "status", "reason", "issue_codes"])
        for a in audit.records:
            reasons = []
            if a.status in clean.exclude_statuses and not (clean.apply_repairs and a.repair_candidate and audit.policy.repair.mode == "apply"):
                reasons.append(f"status={a.status}")
            codes = {i.code for i in a.issues}
            hit = sorted(codes & set(clean.exclude_codes))
            if hit:
                reasons.append("issue=" + ",".join(hit))
            if clean.exclude_min_severity and SEVERITY_RANK[a.max_severity] >= SEVERITY_RANK[clean.exclude_min_severity]:
                reasons.append(f"severity>={clean.exclude_min_severity}")
            if reasons:
                n_x += 1
                wx.writerow([a.record_id, _flat(a.source_id), a.position, a.status, "; ".join(reasons), ";".join(sorted(codes))])
                continue
            if a.status == "ok":
                structure = a.standardized_smiles if clean.use_standardized_structure else a.parsed_smiles
                src = "standardized" if clean.use_standardized_structure else "parsed"
            else:
                structure, src = a.repair_candidate, "repair_candidate"
            n_w += 1
            wo.writerow([a.record_id, _flat(a.source_id), structure, _flat(a.input_structure), src, _flat(a.group_id)] + [_flat(a.fields.get(c)) for c in orig_cols])
    return n_w, n_x


# ------------------------------------------------------------------------------------ HTML

def _svg(smiles: Optional[str], highlight: Optional[List[int]] = None, size=(220, 150)) -> str:
    if not smiles or Chem is None:
        return "<span class=muted>no structure</span>"
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return "<span class=muted>unparseable</span>"
    try:
        d = rdMolDraw2D.MolDraw2DSVG(*size)
        d.DrawMolecule(m, highlightAtoms=[i for i in (highlight or []) if i < m.GetNumAtoms()])
        d.FinishDrawing()
        return d.GetDrawingText().replace("<?xml version='1.0' encoding='iso-8859-1'?>\n", "")
    except Exception:
        return "<span class=muted>depiction failed</span>"


def _hist_png(values: List[float], title: str) -> Optional[str]:
    if not values:
        return None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover
        return None
    fig, ax = plt.subplots(figsize=(3.2, 2.0), dpi=110)
    ax.hist(values, bins=min(30, max(5, len(values) // 5)), color="#4a6fa5")
    ax.set_title(f"{title} (n={len(values)})", fontsize=9); ax.tick_params(labelsize=7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    buf = io.BytesIO(); fig.tight_layout(); fig.savefig(buf, format="png"); plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _e(x) -> str:
    return html.escape("" if x is None else str(x))


def render_html(audit: DatasetAudit, gate: Optional[GateResult] = None, max_rows: int = 500, max_depictions: int = 60) -> str:
    s = audit.summary
    parts: List[str] = []
    parts.append(f"""<!doctype html><html><head><meta charset="utf-8"><title>ChemLitmus audit — {_e(Path(audit.source_file).name)}</title>
<style>body{{font:14px/1.45 system-ui,sans-serif;margin:24px;color:#222}} h1,h2{{margin:18px 0 8px}} table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{border:1px solid #ddd;padding:4px 6px;vertical-align:top;text-align:left}} th{{background:#f3f4f6;position:sticky;top:0}}
.sev-error{{color:#b00020;font-weight:600}} .sev-warning{{color:#9a6700;font-weight:600}} .sev-info{{color:#555}} .muted{{color:#888}}
.grid{{display:flex;flex-wrap:wrap;gap:10px}} .card{{border:1px solid #ddd;border-radius:6px;padding:8px;min-width:240px}} .pass{{color:#0a7f3f;font-weight:700}} .fail{{color:#b00020;font-weight:700}}
input.filter{{width:100%;padding:6px;margin:6px 0;border:1px solid #bbb;border-radius:4px}} svg{{max-width:100%}} code{{background:#f3f4f6;padding:1px 3px}}</style></head><body>""")
    parts.append(f"<h1>ChemLitmus dataset audit</h1><p><b>Input</b> <code>{_e(Path(audit.source_file).name)}</code> · <b>policy</b> {_e(audit.policy.name)} <code>{audit.policy_hash[:16]}…</code> · <b>identity level</b> {_e(s.identity_level)} · <b>schema</b> v{OUTPUT_SCHEMA_VERSION}</p>")
    if gate is not None:
        cls = "pass" if gate.passed else "fail"
        parts.append(f"<p class={cls}>Gate: {'PASS' if gate.passed else 'FAIL'}" + (" — " + "; ".join(_e(v) for v in gate.violations) if gate.violations else "") + "</p>")
    if not s.processing_complete:
        parts.append("<p class=sev-warning>Processing incomplete: some computations failed (see STD_FAILED / unevaluated alerts). Counts below exclude failed items; they are not negatives.</p>")
    # accounting
    parts.append("<h2>Input accounting</h2><table><tr><th>Total</th><th>ok</th><th>empty</th><th>invalid</th><th>unsupported</th><th>error</th><th>annotated</th></tr>"
                 f"<tr><td>{s.n_total}</td><td>{s.n_ok}</td><td>{s.n_empty}</td><td>{s.n_invalid}</td><td>{s.n_unsupported}</td><td>{s.n_error}</td><td>{s.n_annotated}</td></tr></table>")
    if s.dataset_warnings:
        parts.append("<p><b>Dataset-level notes:</b> " + "; ".join(_e(w) for w in s.dataset_warnings) + "</p>")
    # issues summary
    parts.append("<h2>Issues</h2><div class=grid>")
    for code, n in sorted(s.issues_by_code.items(), key=lambda kv: (-SEVERITY_RANK[ISSUE_CATALOGUE.get(kv[0], ('info',))[0]], -kv[1])):
        sev, desc, action = ISSUE_CATALOGUE.get(code, ("info", "", ""))
        parts.append(f"<div class=card><b class=sev-{sev}>{_e(code)}</b> · {n} record(s)<br><span class=muted>{_e(desc)}</span><br><i>{_e(action)}</i></div>")
    parts.append("</div>")
    # identity
    parts.append(f"<h2>Identity</h2><p>{s.n_groups} distinct compounds at level <b>{_e(s.identity_level)}</b>; {s.n_collapsed} records collapse onto an earlier one. Distinct keys per level (not additive; levels nest): "
                 + ", ".join(f"{_e(k)} {v}" for k, v in s.n_groups_by_level.items()) + "</p>")
    multi = [g for g in audit.groups if g.size > 1][:50]
    if multi:
        parts.append("<table><tr><th>group</th><th>size</th><th>records</th><th>differs by</th><th>splits</th></tr>")
        for g in multi:
            parts.append(f"<tr><td>{_e(g.group_id)}</td><td>{g.size}</td><td>{_e(', '.join(x or '' for x in g.source_ids) or ', '.join(g.record_ids))}</td><td>{_e(', '.join(g.differs_by))}</td><td>{_e(', '.join(g.splits))}</td></tr>")
        parts.append("</table>")
    # descriptors
    if s.descriptors:
        parts.append("<h2>Descriptor distributions</h2><div class=grid>")
        for name, ds in s.descriptors.items():
            vals = [a.descriptors[name] for a in audit.records if name in a.descriptors]
            png = _hist_png(vals, name)
            img = f"<img src='data:image/png;base64,{png}'>" if png else ""
            parts.append(f"<div class=card>{img}<br><span class=muted>n={ds.n} min={ds.min:g} median={ds.median:g} max={ds.max:g}</span></div>")
        parts.append("</div>")
    if s.alerts_by_set:
        parts.append("<h2>Structural alerts</h2><p>Records matching at least one alert, by set (denominator {}): ".format(s.alert_denominator)
                     + ", ".join(f"{_e(k)} {v} ({v / max(s.alert_denominator, 1):.1%})" for k, v in s.alerts_by_set.items())
                     + ". Alerts flag assay-interference motifs; they are not evidence of toxicity, activity or synthesisability.</p>")
    # leakage
    if audit.leakage:
        L = audit.leakage
        parts.append(f"<h2>Split leakage</h2><p>Identity level {_e(L['identity_level'])}; records per split: " + ", ".join(f"{_e(k)} {v}" for k, v in L["n_by_split"].items())
                     + f"; unassigned {L['n_unassigned']}. Within-split duplicates: " + ", ".join(f"{_e(k)} {v}" for k, v in L["within_split_duplicates"].items()) + f".<br><span class=muted>{_e(L['note'])}</span></p>")
        for p in L["pairs"]:
            parts.append(f"<h3>{_e(p['reference_split'])} → {_e(p['evaluation_split'])} ({p['n_reference']} vs {p['n_evaluation']})</h3><table><tr><th>class</th><th>eval records</th><th>shared keys</th><th>fraction</th></tr>")
            for o in p["overlap"]:
                parts.append(f"<tr><td>{_e(o['level'])}</td><td>{o['n_eval_records']}</td><td>{o['n_eval_groups']}</td><td>{o['fraction_of_eval']:.1%}</td></tr>")
            parts.append(f"<tr><td class=muted>formula match (not identity)</td><td>{p['formula_matches']}</td><td></td><td></td></tr></table>")
            parts.append(f"<p>Scaffold overlap: {p['scaffold_overlap_records']} of {p['n_evaluation'] - p['n_acyclic_evaluation']} cyclic evaluation records ({p['scaffold_overlap_fraction']:.1%}); {p['n_acyclic_evaluation']} acyclic excluded. "
                         f"Nearest-neighbour {_e(p['neighbour_fingerprint'])} Tanimoto: {_e(p['nearest_neighbour_similarity'])}; {p['n_related_by_similarity']} at or above {p['neighbour_threshold']}. Temporal: {p['temporal_violations']} violation(s) — {_e(p['temporal_note'])}</p>")
    if audit.label_conflicts:
        C = audit.label_conflicts
        parts.append(f"<h2>Label conflicts</h2><p>Endpoint <code>{_e(C['endpoint_field'])}</code> ({_e(C['kind'])}), context {_e(C['context_fields'] or 'none')}, identity level {_e(C['identity_level'])}, tolerance {C['tolerance']}. "
                     f"{C['n_groups_compared']} groups compared, <b>{C['n_conflicts']} in conflict</b> ({C['n_records_in_conflict']} records); censored {C['n_censored']}, missing units {C['n_missing_units']}, unparseable {C['n_unparseable']}.<br><span class=muted>{_e(C['note'])}</span></p>")
        if C["groups"]:
            parts.append("<table><tr><th>group</th><th>records</th><th>spread / labels</th><th>scale</th><th>measurements</th></tr>")
            for g in C["groups"][:50]:
                ms = "; ".join(f"{_e(m['source_id'] or m['record_id'])}: {_e(m['relation'])}{_e(m['raw_value'])} {_e(m['units'] or '')}" for m in g.get("measurements", [])) or _e(g.get("labels"))
                parts.append(f"<tr><td>{_e(g['group_id'])}</td><td>{g['n_records']}</td><td>{_e(g.get('spread') if g.get('spread') is not None else g.get('labels'))}</td><td>{_e(g.get('scale') or '')}</td><td>{ms}</td></tr>")
            parts.append("</table>")
    # standardisation examples
    ex = [a for a in audit.records if a.standardization_changed_identity][:12]
    if ex:
        parts.append("<h2>Standardisation changed identity — examples</h2><table><tr><th>record</th><th>before</th><th>after</th><th>relation</th></tr>")
        for a in ex:
            rel = next((t.identity_relation for t in a.transformations if t.changed_identity), "")
            parts.append(f"<tr><td>{_e(a.source_id or a.record_id)}</td><td>{_svg(a.parsed_smiles)}<br><code>{_e(a.parsed_smiles)}</code></td><td>{_svg(a.standardized_smiles)}<br><code>{_e(a.standardized_smiles)}</code></td><td>{_e(rel)}</td></tr>")
        parts.append("</table>")
    # issue table with depictions
    parts.append("<h2>Issue table</h2><input class=filter placeholder='filter rows (code, id, text)…' oninput=\"var q=this.value.toLowerCase();document.querySelectorAll('#issues tbody tr').forEach(function(r){r.style.display=r.innerText.toLowerCase().indexOf(q)>=0?'':'none'})\">")
    parts.append("<table id=issues><thead><tr><th>severity</th><th>code</th><th>record</th><th>structure</th><th>message</th><th>evidence</th></tr></thead><tbody>")
    by_id = {a.record_id: a for a in audit.records}
    shown = sorted(audit.issues, key=lambda x: (-SEVERITY_RANK[x.severity], x.code))[:max_rows]
    n_dep = 0
    for i in shown:
        a = by_id.get(i.record_id or "")
        dep = ""
        if a is not None and n_dep < max_depictions and a.parsed_smiles:
            hl = i.evidence.get("atoms") if i.code == "ALERT_MATCH" else None
            dep = _svg(a.standardized_smiles or a.parsed_smiles, hl, size=(180, 120)); n_dep += 1
        parts.append(f"<tr><td class=sev-{_e(i.severity)}>{_e(i.severity)}</td><td>{_e(i.code)}</td><td>{_e(i.source_id or i.record_id)}</td><td>{dep}</td><td>{_e(i.message)}</td><td><code>{_e(json.dumps(i.evidence, default=str)[:300])}</code></td></tr>")
    parts.append("</tbody></table>")
    if len(audit.issues) > max_rows:
        parts.append(f"<p class=muted>{len(audit.issues) - max_rows} more issues in issues.csv</p>")
    parts.append("<p class=muted>Descriptor rules and structural alerts do not establish toxicity, ADMET behaviour, biological activity or synthesisability. All counts refer to the denominators stated beside them.</p></body></html>")
    return "\n".join(parts)


# ------------------------------------------------------------------------------------ orchestration

def write_audit_outputs(audit: DatasetAudit, out_dir: str | Path, input_path: Optional[str | Path] = None, gates: Optional[GatePolicy] = None,
                        clean: Optional[CleanPolicy] = None, command: str = "chemlitmus audit") -> Dict[str, str]:
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    paths: Dict[str, Path] = {
        "records": out / "records.csv", "issues": out / "issues.csv", "groups": out / "identity_groups.csv",
        "summary": out / "summary.json", "audit": out / "audit.json", "policy": out / "policy.json", "report": out / "report.html", "manifest": out / "manifest.json",
    }
    gate = evaluate_gates(audit, gates) if gates else None
    write_records_csv(audit, paths["records"]); write_issues_csv(audit, paths["issues"]); write_groups_csv(audit, paths["groups"])
    summary = audit.summary.model_dump(); summary.update({"schema_version": OUTPUT_SCHEMA_VERSION, "policy_hash": audit.policy_hash, "roles": audit.roles, "gate": gate.model_dump() if gate else None})
    paths["summary"].write_text(json.dumps(summary, indent=2, default=str))
    full = audit.model_dump(mode="json"); full["output_schema_version"] = OUTPUT_SCHEMA_VERSION; full["gate"] = gate.model_dump() if gate else None
    paths["audit"].write_text(json.dumps(full, default=str))
    audit.policy.save(paths["policy"])
    paths["report"].write_text(render_html(audit, gate), encoding="utf-8")
    if clean is not None:
        paths["clean"] = out / "clean.csv"; paths["exclusions"] = out / "exclusions.csv"
        write_clean(audit, paths["clean"], paths["exclusions"], clean)
    m = new_manifest(command, identity_level=audit.policy.identity_level, preparation=audit.policy.preparation, alert_sets=audit.policy.alert_sets,
                     fingerprint=f"{audit.policy.fingerprint}/{audit.policy.fingerprint_bits}")
    if input_path and Path(input_path).exists():
        m.add_input(input_path)
    m.policy_hash, m.policy_name, m.output_schema_version = audit.policy_hash, audit.policy.name, OUTPUT_SCHEMA_VERSION
    m.rule_catalogue = {"source": "rdkit FilterCatalog", "sets": ",".join(audit.policy.alert_sets)}
    for k, p in paths.items():
        if k != "manifest":
            m.add_output(p)
    m.save(paths["manifest"])
    return {k: str(p) for k, p in paths.items()}


__all__ = ["GatePolicy", "GateResult", "CleanPolicy", "evaluate_gates", "write_audit_outputs", "render_html", "write_clean", "OUTPUT_SCHEMA_VERSION"]
