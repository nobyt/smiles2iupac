"""Phase 228: formaldehyde retained name (IUPAC 2013 P-66.6.3.1).

C=O → "formaldehyde" is a preferred IUPAC name.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C=O",  "formaldehyde"),
    # regression: other aldehydes still work
    ("CC=O", "acetaldehyde"),
    ("CCC=O", "propanal"),
])
def test_phase228_formaldehyde(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
