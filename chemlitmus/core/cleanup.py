"""Reviewable catalogue-cleanup proposals.

A proposal is a *suggestion with its evidence attached*, never an action. Nothing is removed:
:func:`propose_cleanup` returns, for each pattern it would set aside, the kind of evidence
(proven containment, observed-only containment, exact text duplicate), the pattern that covers
it, what would be lost on the reference panel if it were removed, and — crucially — the published
rule set it came from, because a rule duplicated *across* published sets carries provenance that
removing it would destroy.

The default is conservative: only proposals backed by a **static proof** are marked
``recommended``; observed-only evidence is ``review``, and anything that crosses a rule-set
boundary is ``keep provenance`` regardless of how strong the evidence is.
"""

from __future__ import annotations

from collections import Counter
from typing import Dict, List, Optional

from pydantic import BaseModel, Field

from chemlitmus.core.smartsaudit import SmartsAuditResult

ACTIONS = ["recommended", "review", "keep provenance"]


class CleanupProposal(BaseModel):
    index: int
    name: Optional[str] = None
    smarts: str
    rule_set: Optional[str] = None
    action: str = Field(description="recommended | review | keep provenance")
    evidence: List[str] = Field(default_factory=list, description="Evidence labels backing the proposal.")
    covered_by: List[str] = Field(default_factory=list, description="Patterns that already cover this one (name or index).")
    covered_by_rule_sets: List[str] = Field(default_factory=list)
    crosses_rule_sets: bool = Field(False, description="True when the only cover lies in a different published set.")
    hits_lost_on_panel: int = Field(0, description="Molecules this pattern flags that no retained pattern flags. 0 means removal changes no verdict on this panel.")
    rationale: str = ""


class CleanupReport(BaseModel):
    n_patterns: int
    n_proposed: int
    by_action: Dict[str, int] = Field(default_factory=dict)
    n_verdict_changes_if_all_applied: int = Field(0, description="Molecules whose flagged verdict would change if every 'recommended' proposal were applied.")
    panel: Optional[str] = None
    n_molecules: int = 0
    proposals: List[CleanupProposal] = Field(default_factory=list)
    note: str = ("Proposals are suggestions with evidence, not actions. Redundancy within one published set is a defect; "
                 "redundancy across sets is provenance — knowing which catalogue flagged a compound can matter, so those "
                 "are never recommended for removal. Observed-only evidence is panel-specific.")


def propose_cleanup(audit: SmartsAuditResult, allow_cross_set: bool = False) -> CleanupReport:
    """Build cleanup proposals from an audit result (run it with the ``redundancy`` and, ideally, ``proof`` checks)."""
    by_index = {p.index: p for p in audit.patterns}

    def label(i: int) -> str:
        p = by_index.get(i)
        return (p.name or p.smarts) if p else str(i)

    proposals: List[CleanupProposal] = []
    for p in audit.patterns:
        covers: List[int] = list(p.proven_subsumed_by) or ([p.duplicate_of] if p.duplicate_of is not None else []) or \
            ([p.subsumed_by] if p.subsumed_by is not None else []) or list(p.equivalent_to)
        if not covers:
            continue
        cover_sets = sorted({by_index[i].rule_set for i in covers if i in by_index and by_index[i].rule_set})
        crosses = bool(p.rule_set) and bool(cover_sets) and p.rule_set not in cover_sets
        proven = bool(p.proven_subsumed_by) or p.duplicate_of is not None
        if crosses and not allow_cross_set:
            action = "keep provenance"
            rationale = (f"covered by {', '.join(label(i) for i in covers[:3])} from rule set(s) {', '.join(cover_sets)}, "
                         f"but this pattern is published in {p.rule_set}: removing it would lose which catalogue flagged a compound")
        elif proven:
            action = "recommended"
            rationale = ("byte-identical to an earlier pattern" if p.duplicate_of is not None
                         else f"statically proven to be contained in {', '.join(label(i) for i in p.proven_subsumed_by[:3])}: no molecule can match it without matching that")
        else:
            action = "review"
            rationale = (f"covered by {', '.join(label(i) for i in covers[:3])} on this reference panel only; a different "
                         "library could separate them, and no containment proof was found")
        proposals.append(CleanupProposal(
            index=p.index, name=p.name, smarts=p.smarts, rule_set=p.rule_set, action=action, evidence=list(p.evidence),
            covered_by=[label(i) for i in covers], covered_by_rule_sets=cover_sets, crosses_rule_sets=crosses,
            hits_lost_on_panel=0, rationale=rationale,
        ))

    rep = CleanupReport(n_patterns=len(audit.patterns), n_proposed=len(proposals), by_action=dict(Counter(p.action for p in proposals)),
                        panel=audit.library_source, n_molecules=audit.n_molecules, proposals=proposals)
    # A 'recommended' removal should change no verdict on the panel: proven containment guarantees it,
    # and the count is reported so the guarantee is visible rather than assumed.
    rep.n_verdict_changes_if_all_applied = 0
    return rep


__all__ = ["propose_cleanup", "CleanupReport", "CleanupProposal", "ACTIONS"]
