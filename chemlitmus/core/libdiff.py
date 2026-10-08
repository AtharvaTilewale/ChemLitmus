"""Structure-aware comparison of two compound collections.

A textual ``diff`` of two SMILES files is useless: the same compound is written differently by
different toolkits and database releases, and a salt form, tautomer or stereo assignment can
change without the compound changing. This module compares two collections at a chosen identity
level (see :mod:`chemlitmus.core.identity`) and reports what was added, what was removed, what is
unchanged, and — for compounds present in both — exactly how their representation changed.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Sequence

from pydantic import BaseModel, Field

from chemlitmus.core.identity import (
    IDENTITY_LEVELS,
    IdentityKeys,
    compute_identity,
    describe_difference,
)


class DiffEntry(BaseModel):
    """One compound in the comparison."""

    status: str = Field(description="'added' | 'removed' | 'unchanged' | 'changed'")
    key: str = Field(description="Identity key at the comparison level.")
    smiles_a: List[str] = Field(default_factory=list, description="Input SMILES in collection A sharing this key.")
    smiles_b: List[str] = Field(default_factory=list, description="Input SMILES in collection B sharing this key.")
    count_a: int = 0
    count_b: int = 0
    change: Optional[str] = Field(None, description="For 'changed': what differs between the A and B representations.")


class LibraryDiff(BaseModel):
    """Result of comparing two collections at one identity level."""

    level: str
    n_a: int
    n_b: int
    n_valid_a: int
    n_valid_b: int
    n_keys_a: int = Field(description="Distinct compounds in A at this level.")
    n_keys_b: int = Field(description="Distinct compounds in B at this level.")
    n_added: int = 0
    n_removed: int = 0
    n_unchanged: int = 0
    n_changed: int = Field(0, description="Present in both, but with a different exact representation.")
    changes_by_kind: Dict[str, int] = Field(default_factory=dict)
    multiplicity_changes: int = Field(0, description="Keys present in both whose record count differs.")
    entries: List[DiffEntry] = Field(default_factory=list)
    error: Optional[str] = None

    @property
    def jaccard(self) -> float:
        union = self.n_added + self.n_removed + self.n_unchanged + self.n_changed
        return (self.n_unchanged + self.n_changed) / union if union else 1.0


def _bucket(keys: Sequence[IdentityKeys], level: str) -> Dict[str, List[IdentityKeys]]:
    out: Dict[str, List[IdentityKeys]] = defaultdict(list)
    for k in keys:
        if k.is_valid:
            out[k.key(level)].append(k)
    return out


def diff_libraries(
    smiles_a: Sequence[str],
    smiles_b: Sequence[str],
    level: str = "parent",
    keys_a: Optional[Sequence[IdentityKeys]] = None,
    keys_b: Optional[Sequence[IdentityKeys]] = None,
    include_unchanged: bool = False,
) -> LibraryDiff:
    """Compare two SMILES collections at an identity level.

    Args:
        smiles_a: The reference (older) collection.
        smiles_b: The comparison (newer) collection.
        level: One of :data:`IDENTITY_LEVELS`. ``exact`` compares canonical SMILES as given;
               ``parent`` ignores salts and charges; ``nostereo``/``tautomer``/``skeleton``
               progressively ignore stereo and tautomer differences.
        keys_a, keys_b: Precomputed identity keys, same order as the SMILES.
        include_unchanged: Also list unchanged compounds in ``entries``.

    Returns:
        :class:`LibraryDiff`.
    """
    if level not in IDENTITY_LEVELS:
        raise ValueError(f"Unknown identity level {level!r}. Valid: {IDENTITY_LEVELS}")
    ka = list(keys_a) if keys_a is not None else [compute_identity(s) for s in smiles_a]
    kb = list(keys_b) if keys_b is not None else [compute_identity(s) for s in smiles_b]
    if ka and all(not k.is_valid and k.error == "RDKit is not installed." for k in ka):
        return LibraryDiff(level=level, n_a=len(ka), n_b=len(kb), n_valid_a=0, n_valid_b=0, n_keys_a=0, n_keys_b=0,
                           error="RDKit is not installed.")

    ba, bb = _bucket(ka, level), _bucket(kb, level)
    entries: List[DiffEntry] = []
    kinds: Dict[str, int] = defaultdict(int)
    n_added = n_removed = n_unchanged = n_changed = n_mult = 0

    for key in sorted(set(ba) | set(bb)):
        in_a, in_b = ba.get(key, []), bb.get(key, [])
        if in_a and not in_b:
            n_removed += 1
            entries.append(DiffEntry(status="removed", key=key, smiles_a=[k.input_smiles for k in in_a], count_a=len(in_a)))
            continue
        if in_b and not in_a:
            n_added += 1
            entries.append(DiffEntry(status="added", key=key, smiles_b=[k.input_smiles for k in in_b], count_b=len(in_b)))
            continue
        if len(in_a) != len(in_b):
            n_mult += 1
        exact_a, exact_b = {k.exact for k in in_a}, {k.exact for k in in_b}
        if exact_a == exact_b:
            n_unchanged += 1
            if include_unchanged:
                entries.append(DiffEntry(status="unchanged", key=key,
                                         smiles_a=[k.input_smiles for k in in_a], smiles_b=[k.input_smiles for k in in_b],
                                         count_a=len(in_a), count_b=len(in_b)))
            continue
        # Same compound at this level, different exact representation: say how.
        parts = set()
        for x in in_a:
            for y in in_b:
                if x.exact != y.exact:
                    for part in describe_difference(x, y).replace("; also ", ";").split(";"):
                        if part.strip() and part.strip() not in ("identical", "different compounds"):
                            parts.add(part.strip())
        change = "; ".join(sorted(parts)) if parts else "representation"
        kinds[change] += 1
        n_changed += 1
        entries.append(DiffEntry(status="changed", key=key,
                                 smiles_a=[k.input_smiles for k in in_a], smiles_b=[k.input_smiles for k in in_b],
                                 count_a=len(in_a), count_b=len(in_b), change=change))

    order = {"removed": 0, "added": 1, "changed": 2, "unchanged": 3}
    entries.sort(key=lambda e: (order[e.status], e.key))
    return LibraryDiff(
        level=level,
        n_a=len(ka), n_b=len(kb),
        n_valid_a=sum(1 for k in ka if k.is_valid), n_valid_b=sum(1 for k in kb if k.is_valid),
        n_keys_a=len(ba), n_keys_b=len(bb),
        n_added=n_added, n_removed=n_removed, n_unchanged=n_unchanged, n_changed=n_changed,
        changes_by_kind=dict(kinds), multiplicity_changes=n_mult,
        entries=entries,
    )


__all__ = ["DiffEntry", "LibraryDiff", "diff_libraries"]
