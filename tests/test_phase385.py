"""Phase 385: acrylonitrile is not an IUPAC 2013 preferred name (IUPAC 2013 P-66.6.1.1.1).

The preferred IUPAC name for CH2=CH-CN is prop-2-enenitrile.
'Acrylonitrile' is an acceptable alternative but not the PIN.
Acetonitrile (CC#N) IS a retained preferred name and must still work.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acrylonitrile: systematic name is PIN
    ("C=CC#N",            "prop-2-enenitrile"),
    # regression: acetonitrile retained name still works
    ("CC#N",              "acetonitrile"),
    # regression: other nitriles unaffected
    ("CCC#N",             "propanenitrile"),
    ("CCCC#N",            "butanenitrile"),
    ("N#Cc1ccccc1",       "benzonitrile"),
    # cyano substituent on carboxylic acid (Phase 383)
    ("N#CCC(=O)O",        "3-cyanopropanoic acid"),
])
def test_phase385_acrylonitrile_systematic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
