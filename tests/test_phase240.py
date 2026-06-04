"""Phase 240: acetonitrile retained name (IUPAC 2013 P-66.6.1.1.1).

Acetonitrile is a preferred IUPAC name (PIN) per IUPAC 2013.
Acrylonitrile is NOT a PIN; Phase 385 corrected it to prop-2-enenitrile.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CC#N",    "acetonitrile"),
    ("N#CC",    "acetonitrile"),   # same compound different SMILES
    # acrylonitrile: systematic name is PIN (Phase 385)
    ("C=CC#N",  "prop-2-enenitrile"),
    # regression: longer nitriles use systematic names
    ("CCC#N",   "propanenitrile"),
    ("C=CCC#N", "but-3-enenitrile"),
])
def test_phase240_nitrile_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
