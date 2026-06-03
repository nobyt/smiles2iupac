"""Phase 353: Sulfonohydrazide naming (IUPAC 2013).

Compounds of the form R-S(=O)₂-NH-NH₂ are named as {alkane}sulfonohydrazide.
N'-substituted variants use the N'-prefix convention.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Primary sulfonohydrazides
    ("CS(=O)(=O)NN",            "methanesulfonohydrazide"),
    ("CCS(=O)(=O)NN",           "ethanesulfonohydrazide"),
    ("CCCS(=O)(=O)NN",          "propane-1-sulfonohydrazide"),
    ("CCCCS(=O)(=O)NN",         "butane-1-sulfonohydrazide"),
    ("c1ccc(S(=O)(=O)NN)cc1",   "benzenesulfonohydrazide"),
    # N'-substituted
    ("CS(=O)(=O)NNC",           "N'-methylmethanesulfonohydrazide"),
    ("CCS(=O)(=O)NNC",          "N'-methylethanesulfonohydrazide"),
    ("CS(=O)(=O)NN(C)C",        "N',N'-dimethylmethanesulfonohydrazide"),
    # regressions: sulfonamide and acyl hydrazide unchanged
    ("CS(=O)(=O)N",             "methanesulfonamide"),
    ("CS(=O)(=O)NC",            "N-methylmethanesulfonamide"),
    ("CCS(=O)(=O)NC",           "N-methylethanesulfonamide"),
    ("CC(=O)NN",                "ethanohydrazide"),
    ("CCC(=O)NN",               "propanohydrazide"),
])
def test_phase353_sulfonohydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
