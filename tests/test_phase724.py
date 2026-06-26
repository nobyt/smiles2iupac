"""Phase 724: methyl derivatives of polycyclic aromatic hydrocarbons and ketones:
fluorene, acenaphthylene, acenaphthene, pyrene, fluoranthene, azulene,
coronene, chrysene, perylene, triphenylene, fluoren-9-one,
anthracen-9(10H)-one, naphthalen-1(2H)-one, naphthalen-2(1H)-one.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # fluorene (CH at 1,2,3,4,9)
    ("c1ccc2c(c1)Cc1ccccc1-2",      "fluorene"),
    ("Cc1cccc2c1Cc1ccccc1-2",       "1-methylfluorene"),
    ("Cc1ccc2c(c1)Cc1ccccc1-2",     "2-methylfluorene"),
    ("Cc1ccc2c(c1)-c1ccccc1C2",     "3-methylfluorene"),
    ("Cc1cccc2c1-c1ccccc1C2",       "4-methylfluorene"),
    ("CC1c2ccccc2-c2ccccc21",       "9-methylfluorene"),
    # acenaphthylene (CH at 1,3,4,5; C2 symmetry: 1=2, 3=8, 4=7, 5=6)
    ("C1=Cc2cccc3cccc1c23",         "acenaphthylene"),
    ("CC1=Cc2cccc3cccc1c23",        "1-methylacenaphthylene"),
    ("Cc1ccc2cccc3c2c1C=C3",        "3-methylacenaphthylene"),
    ("Cc1cc2c3c(cccc3c1)C=C2",      "4-methylacenaphthylene"),
    ("Cc1ccc2c3c(cccc13)C=C2",      "5-methylacenaphthylene"),
    # acenaphthene (CH at 1,3,4,5; C2 symmetry: 1=2, 3=8, 4=7, 5=6)
    ("c1cc2c3c(cccc3c1)CC2",        "acenaphthene"),
    ("CC1Cc2cccc3cccc1c23",         "1-methylacenaphthene"),
    ("Cc1ccc2cccc3c2c1CC3",         "3-methylacenaphthene"),
    ("Cc1cc2c3c(cccc3c1)CC2",       "4-methylacenaphthene"),
    ("Cc1ccc2c3c(cccc13)CC2",       "5-methylacenaphthene"),
    # pyrene (CH at 1,2,4; D2h symmetry: 1=3=6=8, 2=7, 4=5=9=10)
    ("c1cc2ccc3cccc4ccc(c1)c2c34",  "pyrene"),
    ("Cc1ccc2ccc3cccc4ccc1c2c34",   "1-methylpyrene"),
    ("Cc1cc2ccc3cccc4ccc(c1)c2c34", "2-methylpyrene"),
    ("Cc1cc2cccc3ccc4cccc1c4c32",   "4-methylpyrene"),
    # fluoranthene (CH at 1,2,3,7,8; C2v symmetry: 1=6, 2=5, 3=4, 7=10, 8=9)
    ("c1ccc2c(c1)-c1cccc3cccc-2c13",    "fluoranthene"),
    ("Cc1ccc2cccc3c2c1-c1ccccc1-3",     "1-methylfluoranthene"),
    ("Cc1cc2c3c(cccc3c1)-c1ccccc1-2",   "2-methylfluoranthene"),
    ("Cc1ccc2c3c(cccc13)-c1ccccc1-2",   "3-methylfluoranthene"),
    ("Cc1cccc2c1-c1cccc3cccc-2c13",     "7-methylfluoranthene"),
    ("Cc1ccc2c(c1)-c1cccc3cccc-2c13",   "8-methylfluoranthene"),
    # azulene (CH at 1,2,4,5,6; C2v symmetry: 1=3, 4=8, 5=7)
    ("c1ccc2cccc-2cc1",             "azulene"),
    ("Cc1ccc2cccccc1-2",            "1-methylazulene"),
    ("Cc1cc2cccccc-2c1",            "2-methylazulene"),
    ("Cc1ccccc2cccc1-2",            "4-methylazulene"),
    ("Cc1cccc2cccc-2c1",            "5-methylazulene"),
    ("Cc1ccc2cccc-2cc1",            "6-methylazulene"),
    # coronene (all CH equivalent → only locant 1)
    ("c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61", "coronene"),
    ("Cc1cc2ccc3ccc4ccc5ccc6ccc1c1c6c5c4c3c21", "1-methylcoronene"),
    # chrysene (CH at 1,2,3,4,5,6; D2h symmetry: 1=7,2=8,3=9,4=10,5=11,6=12)
    ("c1ccc2c(c1)ccc1c3ccccc3ccc21",    "chrysene"),
    ("Cc1cccc2c1ccc1c3ccccc3ccc21",     "1-methylchrysene"),
    ("Cc1ccc2c(ccc3c4ccccc4ccc23)c1",   "2-methylchrysene"),
    ("Cc1ccc2ccc3c4ccccc4ccc3c2c1",     "3-methylchrysene"),
    ("Cc1cccc2ccc3c4ccccc4ccc3c12",     "4-methylchrysene"),
    ("Cc1cc2ccccc2c2ccc3ccccc3c12",     "5-methylchrysene"),
    ("Cc1cc2c3ccccc3ccc2c2ccccc12",     "6-methylchrysene"),
    # perylene (CH at 1,2,3; D2h symmetry: 1=6=7=12, 2=5=8=11, 3=4=9=10)
    ("c1cc2cccc3c4cccc5cccc(c(c1)c23)c54",  "perylene"),
    ("Cc1ccc2cccc3c4cccc5cccc(c1c23)c54",   "1-methylperylene"),
    ("Cc1cc2cccc3c4cccc5cccc(c(c1)c23)c54", "2-methylperylene"),
    ("Cc1ccc2c3cccc4cccc(c5cccc1c52)c43",   "3-methylperylene"),
    # triphenylene (CH at 1,2; D3h symmetry: 1=5=9, 2=6=10)
    ("c1ccc2c(c1)c1ccccc1c1ccccc21",    "triphenylene"),
    ("Cc1cccc2c3ccccc3c3ccccc3c12",     "1-methyltriphenylene"),
    ("Cc1ccc2c3ccccc3c3ccccc3c2c1",     "2-methyltriphenylene"),
    # fluoren-9-one (CH at 1,2,3,4; C2v symmetry: 1=8, 2=7, 3=6, 4=5)
    ("O=C1c2ccccc2-c2ccccc21",          "fluoren-9-one"),
    ("Cc1cccc2c1C(=O)c1ccccc1-2",       "1-methylfluoren-9-one"),
    ("Cc1ccc2c(c1)C(=O)c1ccccc1-2",     "2-methylfluoren-9-one"),
    ("Cc1ccc2c(c1)-c1ccccc1C2=O",       "3-methylfluoren-9-one"),
    ("Cc1cccc2c1-c1ccccc1C2=O",         "4-methylfluoren-9-one"),
    # anthracen-9(10H)-one (CH at 1,2,3,4,10; C2v symmetry: 1=8,2=7,3=6,4=5)
    ("O=C1c2ccccc2Cc2ccccc21",              "anthracen-9(10H)-one"),
    ("Cc1cccc2c1C(=O)c1ccccc1C2",           "1-methylanthracen-9(10H)-one"),
    ("Cc1ccc2c(c1)C(=O)c1ccccc1C2",         "2-methylanthracen-9(10H)-one"),
    ("Cc1ccc2c(c1)Cc1ccccc1C2=O",           "3-methylanthracen-9(10H)-one"),
    ("Cc1cccc2c1Cc1ccccc1C2=O",             "4-methylanthracen-9(10H)-one"),
    ("CC1c2ccccc2C(=O)c2ccccc21",           "10-methylanthracen-9(10H)-one"),
    # naphthalen-1(2H)-one (CH at 2,3,4,5,6,7,8)
    ("O=C1CC=Cc2ccccc21",       "naphthalen-1(2H)-one"),
    ("CC1C=Cc2ccccc2C1=O",      "2-methylnaphthalen-1(2H)-one"),
    ("CC1=Cc2ccccc2C(=O)C1",    "3-methylnaphthalen-1(2H)-one"),
    ("CC1=CCC(=O)c2ccccc21",    "4-methylnaphthalen-1(2H)-one"),
    ("Cc1cccc2c1C=CCC2=O",      "5-methylnaphthalen-1(2H)-one"),
    ("Cc1ccc2c(c1)C=CCC2=O",    "6-methylnaphthalen-1(2H)-one"),
    ("Cc1ccc2c(c1)C(=O)CC=C2",  "7-methylnaphthalen-1(2H)-one"),
    ("Cc1cccc2c1C(=O)CC=C2",    "8-methylnaphthalen-1(2H)-one"),
    # naphthalen-2(1H)-one (CH at 1,3,4,5,6,7,8)
    ("O=C1C=Cc2ccccc2C1",       "naphthalen-2(1H)-one"),
    ("CC1C(=O)C=Cc2ccccc21",    "1-methylnaphthalen-2(1H)-one"),
    ("CC1=Cc2ccccc2CC1=O",      "3-methylnaphthalen-2(1H)-one"),
    ("CC1=CC(=O)Cc2ccccc21",    "4-methylnaphthalen-2(1H)-one"),
    ("Cc1cccc2c1C=CC(=O)C2",    "5-methylnaphthalen-2(1H)-one"),
    ("Cc1ccc2c(c1)C=CC(=O)C2",  "6-methylnaphthalen-2(1H)-one"),
    ("Cc1ccc2c(c1)CC(=O)C=C2",  "7-methylnaphthalen-2(1H)-one"),
    ("Cc1cccc2c1CC(=O)C=C2",    "8-methylnaphthalen-2(1H)-one"),
])
def test_phase724(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
