"""Endpoint and duplicate-label conflicts within identity groups.

Records of the same compound (at the policy identity level, and within the same endpoint
context — target / assay / source when those columns exist) are compared:

* **classification** labels: any disagreement is a conflict;
* **quantitative** measurements: exact (``=``) values are converted to a common unit with the
  conversion recorded; censored values (``<``, ``>``, ``<=``, ``>=``, ``~``, ranges) are kept as
  bounds and never averaged; a conflict is a spread of exact values above the tolerance (in log10
  units for molar quantities, absolute otherwise).

Nothing is resolved automatically: the default action is review.
"""

from __future__ import annotations

import math
import re
from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence

from pydantic import BaseModel, Field

from chemlitmus.core.policy import ChemicalPolicy

# unit -> factor to nM (molar family); other families are compared only within identical units
MOLAR_TO_NM: Dict[str, float] = {"pm": 1e-3, "nm": 1.0, "um": 1e3, "µm": 1e3, "μm": 1e3, "mm": 1e6, "m": 1e9}
CENSOR_RE = re.compile(r"^\s*(<=|>=|<|>|~|=)?\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*(?:-\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?))?\s*$")


class Measurement(BaseModel):
    record_id: str
    source_id: Optional[str] = None
    raw_value: str
    relation: str = Field("=", description="=, <, >, <=, >=, ~ or 'range'")
    value: Optional[float] = Field(None, description="Numeric value in the original units (lower bound for ranges).")
    upper: Optional[float] = None
    units: Optional[str] = None
    value_nm: Optional[float] = Field(None, description="Converted to nM when the unit is molar.")
    log_value: Optional[float] = Field(None, description="-log10(M) when molar (pX scale).")
    conversion: Optional[str] = None
    censored: bool = False
    comparable: bool = Field(True, description="False when units are missing or not convertible to the group's common unit.")
    context: Dict[str, str] = Field(default_factory=dict)


class ConflictGroup(BaseModel):
    group_id: str
    identity_level: str
    key: str
    context: Dict[str, str] = Field(default_factory=dict, description="target / assay / source values shared by the group.")
    kind: str = Field(description="classification | quantitative")
    n_records: int
    labels: Optional[Dict[str, int]] = Field(None, description="Classification: label -> count.")
    n_exact: int = 0
    n_censored: int = 0
    n_incomparable: int = 0
    spread: Optional[float] = Field(None, description="max - min of exact values on the comparison scale.")
    scale: Optional[str] = Field(None, description="'log10 molar' or 'absolute (<units>)'")
    tolerance: Optional[float] = None
    conflict: bool
    measurements: List[Measurement] = Field(default_factory=list)
    suggested_action: str = "review"


class LabelConflictReport(BaseModel):
    endpoint_field: str
    context_fields: List[str]
    identity_level: str
    kind: str
    tolerance: float
    n_groups_compared: int = Field(description="Identity groups with >= 2 records in the same context.")
    n_conflicts: int
    n_records_in_conflict: int
    n_censored: int = 0
    n_missing_units: int = 0
    n_unparseable: int = 0
    groups: List[ConflictGroup] = Field(default_factory=list, description="Conflicting groups only.")
    record_issues: Dict[str, List[Dict[str, Any]]] = Field(default_factory=dict)
    note: str = "Censored values are bounds, not observations, and are never averaged. Groups are compared only within identical endpoint context."


def parse_measurement(raw: Any, units: Optional[str], relation: Optional[str] = None) -> Measurement:
    s = "" if raw is None else str(raw).strip()
    m = Measurement(record_id="", raw_value=s, units=(units or "").strip() or None)
    rel_col = (relation or "").strip()
    mt = CENSOR_RE.match(s)
    if not mt:
        m.comparable = False
        return m
    rel = rel_col or mt.group(1) or "="
    rel = {"'='": "=", "=": "="}.get(rel, rel)
    m.value = float(mt.group(2))
    if mt.group(3) is not None:
        m.upper = float(mt.group(3)); rel = "range"
    m.relation = rel
    m.censored = rel != "="
    if m.units is None:
        m.comparable = False
        return m
    u = m.units.lower()
    if u in MOLAR_TO_NM:
        m.value_nm = m.value * MOLAR_TO_NM[u]
        if m.value_nm > 0:
            m.log_value = round(-math.log10(m.value_nm * 1e-9), 4)
        m.conversion = f"{m.value} {m.units} -> {m.value_nm:g} nM -> p = {m.log_value}"
    return m


def label_conflicts(
    records: Sequence[Dict[str, Any]],
    policy: Optional[ChemicalPolicy] = None,
    endpoint_field: str = "endpoint",
    units_field: Optional[str] = None,
    relation_field: Optional[str] = None,
    context_fields: Sequence[str] = (),
    tolerance: float = 1.0,
    kind: str = "auto",
) -> LabelConflictReport:
    """Detect conflicts.

    Args:
        records: dicts with ``record_id``, ``source_id``, ``key`` (identity key at the policy level) and the raw field values.
        tolerance: log10 units for molar quantities (1.0 = ten-fold), absolute otherwise.
        kind: ``classification`` | ``quantitative`` | ``auto`` (quantitative when every non-empty value parses as a number).
    """
    pol = policy or ChemicalPolicy.preset("parent")
    vals = [str(r.get(endpoint_field, "")).strip() for r in records if str(r.get(endpoint_field, "")).strip()]
    if kind == "auto":
        kind = "quantitative" if vals and all(CENSOR_RE.match(v) for v in vals) else "classification"
    groups: Dict[tuple, List[Dict[str, Any]]] = defaultdict(list)
    issues: Dict[str, List[Dict[str, Any]]] = {}
    n_missing = 0
    if kind == "quantitative" and units_field:
        # units are required for every measurement, grouped or not
        for r in records:
            v = str(r.get(endpoint_field, "")).strip()
            if v and CENSOR_RE.match(v) and not str(r.get(units_field, "")).strip():
                n_missing += 1
                issues.setdefault(r["record_id"], []).append({"code": "UNITS_MISSING", "value": v})
    for r in records:
        if not r.get("key") or not str(r.get(endpoint_field, "")).strip():
            continue
        ctx = tuple((c, str(r.get(c, "")).strip()) for c in context_fields)
        groups[(r["key"], ctx)].append(r)
    out: List[ConflictGroup] = []
    n_compared = n_cens = n_unp = 0
    for n, ((key, ctx), members) in enumerate(groups.items()):
        if len(members) < 2:
            continue
        n_compared += 1
        gid = f"c{n:06d}"
        context = dict(ctx)
        if kind == "classification":
            labels = defaultdict(int)
            for r in members:
                labels[str(r[endpoint_field]).strip()] += 1
            conflict = len(labels) > 1
            cg = ConflictGroup(group_id=gid, identity_level=pol.identity_level, key=key, context=context, kind=kind, n_records=len(members), labels=dict(labels), conflict=conflict)
        else:
            ms: List[Measurement] = []
            for r in members:
                m = parse_measurement(r[endpoint_field], r.get(units_field) if units_field else None, r.get(relation_field) if relation_field else None)
                m.record_id, m.source_id, m.context = r["record_id"], r.get("source_id"), context
                ms.append(m)
                if m.value is None:
                    n_unp += 1
                if m.censored:
                    n_cens += 1
                    issues.setdefault(m.record_id, []).append({"code": "RELATION_CENSORED", "relation": m.relation, "value": m.raw_value})
            exact = [m for m in ms if m.comparable and not m.censored and m.value is not None]
            molar = [m for m in exact if m.log_value is not None]
            if molar and len(molar) == len(exact):
                scale, pts = "log10 molar", [m.log_value for m in molar]
            else:
                units = {m.units.lower() for m in exact if m.units}
                if len(units) == 1:
                    scale, pts = f"absolute ({units.pop()})", [m.value for m in exact]
                else:
                    scale, pts = None, []
                    for m in exact:
                        m.comparable = False
            spread = round(max(pts) - min(pts), 4) if len(pts) >= 2 else None
            conflict = spread is not None and spread > tolerance
            cg = ConflictGroup(group_id=gid, identity_level=pol.identity_level, key=key, context=context, kind=kind, n_records=len(members),
                               n_exact=len(exact), n_censored=sum(1 for m in ms if m.censored), n_incomparable=sum(1 for m in ms if not m.comparable),
                               spread=spread, scale=scale, tolerance=tolerance, conflict=conflict, measurements=ms)
        if cg.conflict:
            out.append(cg)
            for r in members:
                issues.setdefault(r["record_id"], []).append({"code": "LABEL_CONFLICT", "group_id": gid, "kind": kind,
                                                               **({"labels": cg.labels} if cg.labels else {"spread": cg.spread, "scale": cg.scale})})
    return LabelConflictReport(endpoint_field=endpoint_field, context_fields=list(context_fields), identity_level=pol.identity_level, kind=kind, tolerance=tolerance,
                               n_groups_compared=n_compared, n_conflicts=len(out), n_records_in_conflict=sum(g.n_records for g in out),
                               n_censored=n_cens, n_missing_units=n_missing, n_unparseable=n_unp, groups=out, record_issues=issues)


def conflicts_from_annotations(anns, policy: ChemicalPolicy, roles: Dict[str, str], tolerance: float = 1.0) -> LabelConflictReport:
    ctx = [roles[c] for c in ("target", "source") if c in roles]
    recs = []
    for a in anns:
        if not a.identity or not a.identity.is_valid:
            continue
        d = dict(a.fields)
        d.update({"record_id": a.record_id, "source_id": a.source_id, "key": a.identity.key(policy.identity_level)})
        recs.append(d)
    return label_conflicts(recs, policy, endpoint_field=roles["endpoint"], units_field=roles.get("units"), relation_field=roles.get("relation"),
                           context_fields=ctx, tolerance=tolerance)


__all__ = ["label_conflicts", "conflicts_from_annotations", "parse_measurement", "LabelConflictReport", "ConflictGroup", "Measurement", "MOLAR_TO_NM"]
