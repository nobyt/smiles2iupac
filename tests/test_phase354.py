"""Phase 354: Sulfinylhydrazide naming (IUPAC 2013).

Compounds of the form R-S(=O)-NH-NH₂ are named as {alkane}sulfinylhydrazide.
N'-substituted variants use the N'-prefix convention.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Primary sulfinylhydrazides
    ("CS(=O)NN",             "methanesulfinylhydrazide"),
    ("CCS(=O)NN",            "ethanesulfinylhydrazide"),
    ("CCCS(=O)NN",           "propane-1-sulfinylhydrazide"),
    ("c1ccc(S(=O)NN)cc1",    "benzenesulfinylhydrazide"),
    # N'-substituted
    ("CS(=O)NNC",            "N'-methylmethanesulfinylhydrazide"),
    ("CCS(=O)NNC",           "N'-methylethanesulfinylhydrazide"),
    ("CS(=O)NN(C)C",         "N',N'-dimethylmethanesulfinylhydrazide"),
    # regressions: sulfinamide, sulfonohydrazide, sulfonamide unchanged
    ("CS(=O)N",              "methanesulfinamide"),
    ("CS(=O)NC",             "N-methylmethanesulfinamide"),
    ("CS(=O)(=O)NN",         "methanesulfonohydrazide"),
    ("CS(=O)(=O)N",          "methanesulfonamide"),
    ("CC(=O)NN",             "ethanohydrazide"),
])
def test_phase354_sulfinylhydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
