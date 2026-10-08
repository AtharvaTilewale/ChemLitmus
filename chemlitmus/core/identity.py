"""Layered molecular identity.

Two SMILES strings can describe the same compound at one level of resolution and different
compounds at another: a sodium salt and its free acid, two tautomers, a racemate and one of its
enantiomers. This module computes a nested set of identity keys for a molecule and uses them to
group a collection by the level of identity you care about.

Levels, from strictest to loosest (each is a coarsening of the one above it):

================  ================================================================
``exact``         canonical isomeric SMILES of the input, salts and charges included
``parent``        largest organic fragment, neutralised, stereo retained
``tautomer``      parent, tautomer-insensitive (stereo retained)
``nostereo``      parent, stereochemistry removed (tautomer retained)
``skeleton``      parent, stereo removed and tautomer-insensitive
``formula``       molecular formula of the parent
================  ================================================================

``tautomer`` and ``nostereo`` are siblings: both coarsen ``parent`` and both are coarsened by
``skeleton``. Keys for ``tautomer``, ``nostereo`` and ``skeleton`` come from RDKit's
``RegistrationHash`` layers, the same machinery used by compound-registration systems.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from pydantic import BaseModel, Field

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit.Chem import RegistrationHash
    from rdkit.Chem.RegistrationHash import HashLayer

    RDLogger.DisableLog("rdApp.*")
    _RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover
    _RDKIT_AVAILABLE = False


IDENTITY_LEVELS: List[str] = ["exact", "parent", "tautomer", "nostereo", "skeleton", "formula"]
"""Identity levels, strictest first."""

# For a pair that shares a key at the given level but not at ``exact``, the strictest level at
# which they agree tells you what differs between them.
_LEVEL_DIFFERENCE: Dict[str, str] = {
    "parent": "salt, counter-ion or charge form",
    "tautomer": "tautomer",
    "nostereo": "stereochemistry",
    "skeleton": "stereochemistry and tautomer",
    "formula": "constitution (same formula only)",
}


class IdentityKeys(BaseModel):
    """Identity keys for one molecule at every level."""

    input_smiles: str
    is_valid: bool = True
    exact: Optional[str] = None
    parent: Optional[str] = None
    tautomer: Optional[str] = None
    nostereo: Optional[str] = None
    skeleton: Optional[str] = None
    formula: Optional[str] = None
    n_fragments: int = 0
    had_charge: bool = False
    has_stereo: bool = False
    error: Optional[str] = None

    def key(self, level: str) -> Optional[str]:
        if level not in IDENTITY_LEVELS:
            raise ValueError(f"Unknown identity level {level!r}. Valid: {IDENTITY_LEVELS}")
        return getattr(self, level)


class IdentityGroup(BaseModel):
    """A set of input records that share an identity key at the chosen level."""

    level: str
    key: str
    size: int
    indices: List[int] = Field(description="Positions of the member records in the input.")
    smiles: List[str] = Field(description="Input SMILES of the members.")
    distinct_exact: int = Field(description="How many distinct exact-level keys the members have.")
    differs_by: List[str] = Field(
        default_factory=list,
        description="What varies among the members: salt/charge form, tautomer, stereochemistry.",
    )


class IdentityReport(BaseModel):
    """Result of grouping a collection at one identity level."""

    level: str
    n_records: int
    n_valid: int
    n_groups: int = Field(description="Distinct keys at this level (i.e. distinct compounds at this resolution).")
    n_collapsed: int = Field(description="Records that share a key with an earlier record (removable as duplicates at this level).")
    n_groups_by_level: Dict[str, int] = Field(default_factory=dict, description="Distinct keys at every level, for context.")
    groups: List[IdentityGroup] = Field(default_factory=list, description="Groups with more than one member, largest first.")
    keys: List[IdentityKeys] = Field(default_factory=list)
    error: Optional[str] = None


def _parent_mol(mol: "Chem.Mol") -> "Chem.Mol":
    chooser = rdMolStandardize.LargestFragmentChooser()
    uncharger = rdMolStandardize.Uncharger()
    p = chooser.choose(mol)
    p = uncharger.uncharge(p)
    return p


def compute_identity(smiles: str) -> IdentityKeys:
    """Compute identity keys at every level for one SMILES string."""
    if not _RDKIT_AVAILABLE:
        return IdentityKeys(input_smiles=smiles, is_valid=False, error="RDKit is not installed.")
    s = (smiles or "").strip()
    mol = Chem.MolFromSmiles(s) if s else None
    if mol is None:
        return IdentityKeys(input_smiles=smiles, is_valid=False, error="Invalid SMILES: could not be parsed by RDKit.")
    try:
        exact = Chem.MolToSmiles(mol)
        parent = _parent_mol(mol)
        layers = RegistrationHash.GetMolLayers(parent)
        si = Chem.FindPotentialStereo(parent)
        return IdentityKeys(
            input_smiles=smiles,
            exact=exact,
            parent=layers[HashLayer.CANONICAL_SMILES],
            tautomer=layers[HashLayer.TAUTOMER_HASH],
            nostereo=layers[HashLayer.NO_STEREO_SMILES],
            skeleton=layers[HashLayer.NO_STEREO_TAUTOMER_HASH],
            formula=layers[HashLayer.FORMULA],
            n_fragments=len(Chem.GetMolFrags(mol)),
            had_charge=any(a.GetFormalCharge() != 0 for a in mol.GetAtoms()),
            has_stereo=any(e.specified == Chem.StereoSpecified.Specified for e in si),
        )
    except Exception as exc:  # pragma: no cover - defensive
        return IdentityKeys(input_smiles=smiles, is_valid=False, error=f"Identity computation failed: {exc}")


def strictest_shared_level(a: IdentityKeys, b: IdentityKeys) -> Optional[str]:
    """The strictest identity level at which two molecules share a key, or ``None``."""
    for level in IDENTITY_LEVELS:
        ka, kb = a.key(level), b.key(level)
        if ka is not None and ka == kb:
            return level
    return None


def describe_difference(a: IdentityKeys, b: IdentityKeys) -> str:
    """Plain-language description of what differs between two molecules."""
    level = strictest_shared_level(a, b)
    if level is None:
        return "different compounds"
    if level == "exact":
        return "identical"
    desc = _LEVEL_DIFFERENCE[level]
    if level != "parent" and (a.n_fragments != b.n_fragments or a.had_charge != b.had_charge):
        desc += "; also salt, counter-ion or charge form"
    return desc


def _differs_by(members: Sequence[IdentityKeys], level: str) -> List[str]:
    """What varies within a group whose members are identical at ``level``.

    Compares how many distinct keys the members have at each stricter level: if collapsing a
    distinction (stereo, tautomer, salt form) merges keys, that distinction varies in the group.
    """
    n = {lv: len({m.key(lv) for m in members}) for lv in IDENTITY_LEVELS}
    out: List[str] = []
    if n["exact"] > n["parent"]:
        out.append("salt, counter-ion or charge form")
    if n["tautomer"] < n["parent"] or n["skeleton"] < n["nostereo"]:
        out.append("tautomer")
    if n["nostereo"] < n["parent"] or n["skeleton"] < n["tautomer"]:
        out.append("stereochemistry")
    if level == "formula" and n["skeleton"] > 1:
        out.append("constitution")
    return out


def group_by_identity(
    smiles: Iterable[str],
    level: str = "parent",
    keys: Optional[Sequence[IdentityKeys]] = None,
) -> IdentityReport:
    """Group a collection of SMILES by identity at the chosen level.

    Args:
        smiles: Input SMILES strings.
        level: One of :data:`IDENTITY_LEVELS`.
        keys: Precomputed keys (same order as ``smiles``) to avoid recomputation.

    Returns:
        :class:`IdentityReport` with every multi-member group and per-record keys.
    """
    if level not in IDENTITY_LEVELS:
        raise ValueError(f"Unknown identity level {level!r}. Valid: {IDENTITY_LEVELS}")
    smiles_list = list(smiles)
    if not _RDKIT_AVAILABLE:
        return IdentityReport(level=level, n_records=len(smiles_list), n_valid=0, n_groups=0, n_collapsed=0,
                              error="RDKit is not installed.")
    all_keys = list(keys) if keys is not None else [compute_identity(s) for s in smiles_list]
    valid = [(i, k) for i, k in enumerate(all_keys) if k.is_valid]

    buckets: Dict[str, List[int]] = defaultdict(list)
    for i, k in valid:
        buckets[k.key(level)].append(i)

    groups: List[IdentityGroup] = []
    for key, idxs in buckets.items():
        if len(idxs) < 2:
            continue
        members = [all_keys[i] for i in idxs]
        groups.append(IdentityGroup(
            level=level, key=key, size=len(idxs), indices=idxs,
            smiles=[all_keys[i].input_smiles for i in idxs],
            distinct_exact=len({m.exact for m in members}),
            differs_by=_differs_by(members, level),
        ))
    groups.sort(key=lambda g: (-g.size, g.key))

    by_level = {lv: len({k.key(lv) for _, k in valid}) for lv in IDENTITY_LEVELS}
    return IdentityReport(
        level=level,
        n_records=len(smiles_list),
        n_valid=len(valid),
        n_groups=len(buckets),
        n_collapsed=len(valid) - len(buckets),
        n_groups_by_level=by_level,
        groups=groups,
        keys=all_keys,
    )


__all__ = [
    "IDENTITY_LEVELS",
    "IdentityKeys",
    "IdentityGroup",
    "IdentityReport",
    "compute_identity",
    "group_by_identity",
    "strictest_shared_level",
    "describe_difference",
]
