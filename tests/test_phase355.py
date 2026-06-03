"""Phase 355: Sulfonyl azide naming (IUPAC 2013).

Compounds of the form R-S(=O)₂-N₃ are named as {alkane}sulfonyl azide.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Alkyl sulfonyl azides
    ("CS(=O)(=O)N=[N+]=[N-]",              "methanesulfonyl azide"),
    ("CCS(=O)(=O)N=[N+]=[N-]",             "ethanesulfonyl azide"),
    ("CCCS(=O)(=O)N=[N+]=[N-]",            "propane-1-sulfonyl azide"),
    ("CCCCS(=O)(=O)N=[N+]=[N-]",           "butane-1-sulfonyl azide"),
    # Aryl sulfonyl azide
    ("c1ccc(S(=O)(=O)N=[N+]=[N-])cc1",     "benzenesulfonyl azide"),
    # regressions: sulfonamide, sulfonohydrazide, acyl azide unchanged
    ("CS(=O)(=O)N",                         "methanesulfonamide"),
    ("CS(=O)(=O)NN",                        "methanesulfonohydrazide"),
    ("CS(=O)(=O)NC",                        "N-methylmethanesulfonamide"),
    ("CC(=O)N=[N+]=[N-]",                  "acetyl azide"),
])
def test_phase355_sulfonyl_azide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
