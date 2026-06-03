"""Phase 347: O-substituted oxime (oxime ether) naming (IUPAC 2013).

Compounds of the form C=N-O-R are named as O-alkyl{parent}al oxime or
O-alkyl{parent}an-N-one oxime depending on whether the parent carbonyl
is an aldehyde (aldoxime) or ketone (ketoxime) type.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # O-alkyl aldoximes (aldehyde-type)
    ("C=NOC",                    "O-methylmethanal oxime"),
    ("CC=NOC",                   "O-methylethanal oxime"),
    ("CC=NOCC",                  "O-ethylethanal oxime"),
    ("CCCC=NOC",                 "O-methylbutanal oxime"),
    # O-alkyl ketoximes (ketone-type)
    ("CC(=NOC)C",                "O-methylpropan-2-one oxime"),
    ("CCC(=NOC)C",               "O-methylbutan-2-one oxime"),
    ("CC(=NOCC)C",               "O-ethylpropan-2-one oxime"),
    # regressions: regular oximes unchanged
    ("CC=NO",                    "ethanal oxime"),
    ("CC(=NO)C",                 "propan-2-one oxime"),
])
def test_phase347_o_substituted_oxime(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
