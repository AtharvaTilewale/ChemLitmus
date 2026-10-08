"""Tests for chemlitmus.core.augment"""
from chemlitmus.core.augment import augment_smiles

def test_augment_valid():
    # A molecule with multiple canonical SMILES possibilities
    result = augment_smiles("CC(=O)OC1=CC=CC=C1C(=O)O", num_augmentations=5)
    assert result.success is True
    assert len(result.augmented_smiles) > 1

def test_augment_invalid():
    result = augment_smiles("invalid", num_augmentations=5)
    assert result.success is False
    assert result.error is not None
