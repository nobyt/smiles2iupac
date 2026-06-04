"""Phase 419: Naphthalen-1(2H)-one, naphthalen-2(1H)-one,
and 1,3-benzodioxol-2-one.

IUPAC 2013: systematic/retained names for partially reduced naphthalenone
tautomers and the cyclic catechol carbonate.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # naphthalen-1(2H)-one — C1=O, C2 is CH2
    ("O=C1CC=Cc2ccccc21",             "naphthalen-1(2H)-one"),
    # naphthalen-2(1H)-one — C2=O, C1 is CH2
    ("O=C1C=Cc2ccccc2C1",             "naphthalen-2(1H)-one"),
    # 1,3-benzodioxol-2-one — cyclic carbonate of catechol
    ("O=c1oc2ccccc2o1",               "1,3-benzodioxol-2-one"),
    # regression: naphthalene unchanged
    ("c1ccc2ccccc2c1",                 "naphthalene"),
    # regression: 3,4-dihydronaphthalen-1(2H)-one unchanged (Phase 410)
    ("O=C1CCCc2ccccc21",              "3,4-dihydronaphthalen-1(2H)-one"),
    # regression: 1,3-benzodioxole unchanged
    ("c1ccc2c(c1)OCO2",               "1,3-benzodioxole"),
    # regression: coumarin unchanged
    ("O=c1ccc2ccccc2o1",              "coumarin"),
    # regression: benzene unchanged
    ("c1ccccc1",                        "benzene"),
])
def test_phase419_naphthalenones_benzodioxolone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
