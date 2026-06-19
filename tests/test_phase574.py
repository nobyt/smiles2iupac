"""Phase 574: Substituted fluoren-9-one and 1H-phenalene naming.
fluoren-9-one: C2v-symmetric, 14 atoms, 4 unique pairs: {1,8}, {2,7}, {3,6}, {4,5}.
1H-phenalene: no symmetry, 13 atoms, 9 unique positions (1-9).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # fluoren-9-one
    ("O=C1c2ccccc2-c2ccccc21",               "fluoren-9-one"),
    ("Cc1cccc2c1C(=O)c1ccccc1-2",            "1-methylfluoren-9-one"),
    ("Cc1ccc2c(c1)C(=O)c1ccccc1-2",          "2-methylfluoren-9-one"),
    ("Cc1ccc2c(c1)-c1ccccc1C2=O",            "3-methylfluoren-9-one"),
    ("Cc1cccc2c1-c1ccccc1C2=O",              "4-methylfluoren-9-one"),
    # 1H-phenalene
    ("C1=Cc2cccc3cccc(c23)C1",               "1H-phenalene"),
    ("CC1C=Cc2cccc3cccc1c23",                "1-methyl-1H-phenalene"),
    ("CC1=Cc2cccc3cccc(c23)C1",              "2-methyl-1H-phenalene"),
    ("CC1=CCc2cccc3cccc1c23",                "3-methyl-1H-phenalene"),
    ("Cc1ccc2cccc3c2c1C=CC3",               "4-methyl-1H-phenalene"),
    ("Cc1cc2c3c(cccc3c1)CC=C2",             "5-methyl-1H-phenalene"),
    ("Cc1ccc2c3c(cccc13)CC=C2",             "6-methyl-1H-phenalene"),
    ("Cc1ccc2c3c(cccc13)C=CC2",             "7-methyl-1H-phenalene"),
    ("Cc1cc2c3c(cccc3c1)C=CC2",             "8-methyl-1H-phenalene"),
    ("Cc1ccc2cccc3c2c1CC=C3",               "9-methyl-1H-phenalene"),
])
def test_phase574_fluorenone_phenalene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
