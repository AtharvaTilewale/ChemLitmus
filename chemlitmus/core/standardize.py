"""SMILES Standardization Pipeline using RDKit MolStandardize.

All operations are fully offline. No network access required.

Pipeline steps (applied in order):
    1. fragment  — Salt stripping: keep the largest organic fragment.
    2. neutralize — Neutralize charged atoms where chemically reasonable.
    3. tautomer  — Canonicalize tautomers to a single consistent form.
    4. canonical — Generate canonical SMILES string.
"""

from __future__ import annotations

from typing import List, Optional

from pydantic import BaseModel, Field

try:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

# Valid step names (order matters for pipeline execution)
VALID_STEPS: List[str] = ["fragment", "neutralize", "tautomer", "canonical"]


class StepResult(BaseModel):
    """Result for a single standardization step."""

    step: str
    input_smiles: str
    output_smiles: str
    changed: bool = Field(default=False)
    note: str = ""


class StandardizeResult(BaseModel):
    """Complete result of the SMILES standardization pipeline."""

    input_smiles: str
    output_smiles: Optional[str] = None
    steps_applied: List[str] = Field(default_factory=list)
    step_results: List[StepResult] = Field(default_factory=list)
    changed: bool = False
    error: Optional[str] = None

    @property
    def success(self) -> bool:
        return self.error is None and self.output_smiles is not None


def standardize_smiles(
    smiles: str,
    steps: Optional[List[str]] = None,
) -> StandardizeResult:
    """Standardize a SMILES string through a configurable pipeline.

    Args:
        smiles:  Input SMILES string.
        steps:   List of step names to apply. Defaults to all steps in order.
                 Valid values: "fragment", "neutralize", "tautomer", "canonical".
                 Use ["all"] or None to apply every step.

    Returns:
        StandardizeResult with the standardized SMILES and per-step details.

    Raises:
        ValueError: If an unknown step name is provided.
    """
    if not _RDKIT_AVAILABLE:
        return StandardizeResult(
            input_smiles=smiles,
            error="RDKit is not installed. Install it with: pip install rdkit",
        )

    # Resolve steps
    if steps is None or (len(steps) == 1 and steps[0] == "all"):
        steps_to_run = VALID_STEPS[:]
    else:
        invalid = [s for s in steps if s not in VALID_STEPS]
        if invalid:
            raise ValueError(
                f"Unknown step(s): {invalid}. Valid steps: {VALID_STEPS}"
            )
        steps_to_run = [s for s in VALID_STEPS if s in steps]

    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return StandardizeResult(
            input_smiles=smiles,
            error=f"Invalid SMILES: could not parse '{smiles}'",
        )

    current_smiles = Chem.MolToSmiles(mol)
    step_results: List[StepResult] = []

    for step_name in steps_to_run:
        step_in = current_smiles
        try:
            mol, current_smiles = _apply_step(step_name, mol, current_smiles)
        except Exception as exc:
            return StandardizeResult(
                input_smiles=smiles,
                output_smiles=current_smiles,
                steps_applied=[r.step for r in step_results],
                step_results=step_results,
                changed=current_smiles != Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) if smiles else False,
                error=f"Step '{step_name}' failed: {exc}",
            )

        step_results.append(
            StepResult(
                step=step_name,
                input_smiles=step_in,
                output_smiles=current_smiles,
                changed=(current_smiles != step_in),
            )
        )

    original_canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    overall_changed = current_smiles != original_canonical

    return StandardizeResult(
        input_smiles=smiles,
        output_smiles=current_smiles,
        steps_applied=steps_to_run,
        step_results=step_results,
        changed=overall_changed,
    )


def _apply_step(step_name: str, mol, current_smiles: str):
    """Apply a single standardization step. Returns (mol, smiles)."""
    if step_name == "fragment":
        chooser = rdMolStandardize.LargestFragmentChooser()
        mol = chooser.choose(mol)
        smiles = Chem.MolToSmiles(mol)

    elif step_name == "neutralize":
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        smiles = Chem.MolToSmiles(mol)

    elif step_name == "tautomer":
        enumerator = rdMolStandardize.TautomerEnumerator()
        mol = enumerator.Canonicalize(mol)
        smiles = Chem.MolToSmiles(mol)

    elif step_name == "canonical":
        smiles = Chem.MolToSmiles(mol)

    else:
        raise ValueError(f"Unknown step: {step_name}")

    return mol, smiles
