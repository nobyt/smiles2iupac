"""Phase 240: acetonitrile and acrylonitrile retained names (IUPAC 2013 P-66.6.1.1.1).

Both are preferred IUPAC names (PINs) per IUPAC 2013.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CC#N",    "acetonitrile"),
    ("N#CC",    "acetonitrile"),   # same compound different SMILES
    ("C=CC#N",  "acrylonitrile"),
    # regression: longer nitriles use systematic names
    ("CCC#N",   "propanenitrile"),
    ("C=CCC#N", "but-3-enenitrile"),
])
def test_phase240_nitrile_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
