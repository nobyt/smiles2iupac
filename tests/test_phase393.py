"""Phase 393: Peroxy ester naming – {alkyl} {acid}aneperoxoate (IUPAC 2013 P-65.1.5).

CC(=O)OOC was previously undetected as a functional group and fell through to
a generic chain naming, giving '1-oxo-1-oxyethane'.  A new peroxy_ester group
type (priority 95, between peroxyacid=99 and ester=90) detects the pattern
C(=O)-O-O-C and names it correctly.

Acyl peroxides (C(=O)-O-O-C(=O)) are explicitly excluded and still handled
by their own namer.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # basic peroxy esters
    ("CC(=O)OOC",         "methyl ethaneperoxoate"),
    ("CC(=O)OOCC",        "ethyl ethaneperoxoate"),
    ("CCC(=O)OOC",        "methyl propaneperoxoate"),
    ("CCCC(=O)OOCC",      "ethyl butaneperoxoate"),
    # regression: peroxyacid still works
    ("CC(=O)OO",          "ethaneperoxoic acid"),
    ("CCC(=O)OO",         "propaneperoxoic acid"),
    # regression: hydroperoxide still works
    ("COO",               "methyl hydroperoxide"),
    ("CCOO",              "ethyl hydroperoxide"),
    # regression: plain ester still works
    ("CC(=O)OC",          "methyl acetate"),
    ("CC(=O)OCC",         "ethyl acetate"),
    # regression: acyl peroxide still works
    ("CC(=O)OOC(=O)C",    "diethanoyl peroxide"),
])
def test_phase393_peroxyester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
