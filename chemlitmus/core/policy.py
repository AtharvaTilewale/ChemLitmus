"""Versioned chemical policy.

Every chemical decision ChemLitmus makes on a dataset — which fragment is the parent, whether
charges are neutralised, how tautomers and stereo are treated, which preparation SMARTS are
matched under, which fingerprint is used, whether mechanical repairs may be applied — is a
*choice*. :class:`ChemicalPolicy` makes the choice explicit, validates it, and hashes it so a
result can state exactly which policy produced it. Two presets are provided:

* ``conservative`` — keep every component and charge; identity at the ``exact`` level; no repairs.
* ``parent`` — largest organic fragment, neutralised; identity at the ``parent`` level; repairs
  reported but not applied.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict, List, Optional

from pydantic import BaseModel, Field, field_validator

POLICY_SCHEMA_VERSION = "1"
PRESETS = ["conservative", "parent"]


class RepairPolicy(BaseModel):
    mode: str = Field("report", description="'none' (do not attempt), 'report' (attempt, keep original, list candidate), 'apply' (replace the structure with the candidate; original retained in provenance)")
    allowed: List[str] = Field(default_factory=lambda: ["characters", "parentheses", "rings", "brackets", "aromaticity"],
                               description="Problem classes whose mechanical repairs may be considered.")

    @field_validator("mode")
    @classmethod
    def _mode(cls, v):
        if v not in ("none", "report", "apply"):
            raise ValueError("repair.mode must be none, report or apply")
        return v


class SeverityPolicy(BaseModel):
    """Override the default severity of issue codes (info | warning | error)."""

    overrides: Dict[str, str] = Field(default_factory=dict)

    @field_validator("overrides")
    @classmethod
    def _sev(cls, v):
        for k, s in v.items():
            if s not in ("info", "warning", "error", "ignore"):
                raise ValueError(f"severity for {k} must be info, warning, error or ignore")
        return v


class ChemicalPolicy(BaseModel):
    """The resolved set of chemical choices for one run."""

    schema_version: str = POLICY_SCHEMA_VERSION
    name: str = Field("parent", description="Free label; presets are 'conservative' and 'parent'.")
    # representations
    allow_cxsmiles: bool = True
    # fragments / charges
    fragment: str = Field("largest_organic", description="'keep_all' or 'largest_organic' (RDKit LargestFragmentChooser).")
    neutralize: bool = Field(True, description="Apply RDKit Uncharger to the parent.")
    flag_multiple_organic_components: bool = Field(True, description="Flag records whose components include more than one substantial organic fragment (mixtures), for review.")
    organic_component_min_heavy_atoms: int = Field(6, description="A component with at least this many heavy atoms counts as substantial.")
    # tautomer / stereo / isotope
    tautomer: str = Field("keep", description="'keep' or 'canonicalize' (RDKit TautomerEnumerator canonical form) when producing the standardised structure.")
    stereo: str = Field("keep", description="'keep' or 'strip' when producing the standardised structure.")
    isotopes: str = Field("keep", description="'keep' or 'strip'.")
    # identity
    identity_level: str = Field("parent", description="Level used for duplicate / leakage grouping: exact, parent, tautomer, nostereo, skeleton.")
    # matching
    preparation: str = Field("implicit-h", description="Molecule preparation for SMARTS matching: implicit-h, explicit-h, kekule.")
    use_chirality_in_matching: bool = False
    alert_sets: List[str] = Field(default_factory=lambda: ["PAINS"], description="RDKit FilterCatalog sets to apply (PAINS, PAINS_A/B/C, BRENK, NIH, ZINC, CHEMBL_*).")
    # fingerprints
    fingerprint: str = Field("ecfp4", description="ecfp4 | ecfp6 | fcfp4 | rdkit | atompair | torsion | maccs")
    fingerprint_bits: int = 2048
    similarity_threshold: float = Field(0.7, description="Tanimoto at or above which a nearest neighbour is reported as 'related'.")
    # repairs and severities
    repair: RepairPolicy = Field(default_factory=RepairPolicy)
    severity: SeverityPolicy = Field(default_factory=SeverityPolicy)

    @field_validator("fragment")
    @classmethod
    def _frag(cls, v):
        if v not in ("keep_all", "largest_organic"):
            raise ValueError("fragment must be keep_all or largest_organic")
        return v

    @field_validator("tautomer", "stereo", "isotopes")
    @classmethod
    def _keep_or(cls, v, info):
        allowed = {"tautomer": ("keep", "canonicalize"), "stereo": ("keep", "strip"), "isotopes": ("keep", "strip")}[info.field_name]
        if v not in allowed:
            raise ValueError(f"{info.field_name} must be one of {allowed}")
        return v

    @field_validator("identity_level")
    @classmethod
    def _level(cls, v):
        if v not in ("exact", "parent", "tautomer", "nostereo", "skeleton"):
            raise ValueError("identity_level must be exact, parent, tautomer, nostereo or skeleton (formula is not an identity)")
        return v

    @field_validator("preparation")
    @classmethod
    def _prep(cls, v):
        if v not in ("implicit-h", "explicit-h", "kekule"):
            raise ValueError("preparation must be implicit-h, explicit-h or kekule")
        return v

    @field_validator("fingerprint")
    @classmethod
    def _fp(cls, v):
        if v not in ("ecfp4", "ecfp6", "fcfp4", "rdkit", "atompair", "torsion", "maccs"):
            raise ValueError("unknown fingerprint")
        return v

    # ---- hashing / io ----------------------------------------------------------------
    def canonical_json(self) -> str:
        return json.dumps(self.model_dump(mode="json"), sort_keys=True, separators=(",", ":"))

    @property
    def hash(self) -> str:
        """SHA-256 of the canonical JSON; identical policies hash identically."""
        return hashlib.sha256(self.canonical_json().encode()).hexdigest()

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.model_dump(mode="json"), indent=2) + "\n")

    @classmethod
    def load(cls, path: str | Path) -> "ChemicalPolicy":
        return cls.model_validate(json.loads(Path(path).read_text()))

    @classmethod
    def preset(cls, name: str) -> "ChemicalPolicy":
        if name == "conservative":
            return cls(name="conservative", fragment="keep_all", neutralize=False, identity_level="exact",
                       repair=RepairPolicy(mode="none"))
        if name == "parent":
            return cls(name="parent")
        raise ValueError(f"Unknown preset {name!r}. Available: {PRESETS}")


def resolve_policy(config: Optional[str | Path | dict | ChemicalPolicy]) -> ChemicalPolicy:
    """Policy from a preset name, a JSON file, a dict, an instance, or ``None`` (parent preset)."""
    if config is None:
        return ChemicalPolicy.preset("parent")
    if isinstance(config, ChemicalPolicy):
        return config
    if isinstance(config, dict):
        return ChemicalPolicy.model_validate(config)
    s = str(config)
    if s in PRESETS:
        return ChemicalPolicy.preset(s)
    return ChemicalPolicy.load(s)


__all__ = ["ChemicalPolicy", "RepairPolicy", "SeverityPolicy", "resolve_policy", "PRESETS", "POLICY_SCHEMA_VERSION"]
