"""Phase 333: E/Z in sulfonate/sulfinate ester and peroxoic acid; double bond in ester (IUPAC 2013).

Sulfonate/sulfinate esters with unsaturated acid chains now get E/Z prefix.
Peroxoic acids with unsaturated chains now get proper unsaturated name + E/Z.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # sulfonate ester E/Z
    ("C/C=C/CS(=O)(=O)OC",   "(2E)-methyl but-2-ene-1-sulfonate"),
    (r"C/C=C\CS(=O)(=O)OC",  "(2Z)-methyl but-2-ene-1-sulfonate"),
    # sulfinate ester E/Z
    ("C/C=C/CS(=O)OC",       "(2E)-methyl but-2-ene-1-sulfinate"),
    # peroxoic acid with double bond
    ("C/C=C/C(=O)OO",        "(2E)-but-2-eneperoxoic acid"),
    (r"C/C=C\C(=O)OO",       "(2Z)-but-2-eneperoxoic acid"),
    # regressions: saturated sulfonate/sulfinate esters
    ("CS(=O)(=O)OC",         "methyl methanesulfonate"),
    ("CS(=O)OC",             "methyl methanesulfinate"),
    ("CCCS(=O)(=O)OC",       "methyl propane-1-sulfonate"),
    # regressions: saturated peroxoic acids
    ("CC(=O)OO",             "ethaneperoxoic acid"),
    ("CCC(=O)OO",            "propaneperoxoic acid"),
])
def test_phase333_sulfonate_ester_peroxoic_ez(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
