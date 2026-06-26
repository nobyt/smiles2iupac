"""Phase 719: methyl derivatives of quinolone, isoquinolinone, cinnolinone,
phthalazinone, quinazolinone, quinazolinedione, and quinoxalinone series.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # quinolin-2(1H)-one (CH at 3,4,5,6,7,8)
    ("O=c1ccc2ccccc2[nH]1",         "quinolin-2(1H)-one"),
    ("Cc1cc2ccccc2[nH]c1=O",        "3-methylquinolin-2(1H)-one"),
    ("Cc1cc(=O)[nH]c2ccccc12",      "4-methylquinolin-2(1H)-one"),
    ("Cc1cccc2[nH]c(=O)ccc12",      "5-methylquinolin-2(1H)-one"),
    ("Cc1ccc2[nH]c(=O)ccc2c1",      "6-methylquinolin-2(1H)-one"),
    ("Cc1ccc2ccc(=O)[nH]c2c1",      "7-methylquinolin-2(1H)-one"),
    ("Cc1cccc2ccc(=O)[nH]c12",      "8-methylquinolin-2(1H)-one"),
    # quinolin-4(1H)-one (CH at 2,3,5,6,7,8)
    ("O=c1cc[nH]c2ccccc12",         "quinolin-4(1H)-one"),
    ("Cc1cc(=O)c2ccccc2[nH]1",      "2-methylquinolin-4(1H)-one"),
    ("Cc1c[nH]c2ccccc2c1=O",        "3-methylquinolin-4(1H)-one"),
    ("Cc1cccc2[nH]ccc(=O)c12",      "5-methylquinolin-4(1H)-one"),
    ("Cc1ccc2[nH]ccc(=O)c2c1",      "6-methylquinolin-4(1H)-one"),
    ("Cc1ccc2c(=O)cc[nH]c2c1",      "7-methylquinolin-4(1H)-one"),
    ("Cc1cccc2c(=O)cc[nH]c12",      "8-methylquinolin-4(1H)-one"),
    # isoquinolin-1(2H)-one (CH at 3,4,5,6,7,8)
    ("O=c1[nH]ccc2ccccc12",         "isoquinolin-1(2H)-one"),
    ("Cc1cc2ccccc2c(=O)[nH]1",      "3-methylisoquinolin-1(2H)-one"),
    ("Cc1c[nH]c(=O)c2ccccc12",      "4-methylisoquinolin-1(2H)-one"),
    ("Cc1cccc2c(=O)[nH]ccc12",      "5-methylisoquinolin-1(2H)-one"),
    ("Cc1ccc2c(=O)[nH]ccc2c1",      "6-methylisoquinolin-1(2H)-one"),
    ("Cc1ccc2cc[nH]c(=O)c2c1",      "7-methylisoquinolin-1(2H)-one"),
    ("Cc1cccc2cc[nH]c(=O)c12",      "8-methylisoquinolin-1(2H)-one"),
    # isoquinolin-3(2H)-one (CH at 1,4,5,6,7,8)
    ("O=c1cc2ccccc2c[nH]1",         "isoquinolin-3(2H)-one"),
    ("Cc1[nH]c(=O)cc2ccccc12",      "1-methylisoquinolin-3(2H)-one"),
    ("Cc1c(=O)[nH]cc2ccccc12",      "4-methylisoquinolin-3(2H)-one"),
    ("Cc1cccc2c[nH]c(=O)cc12",      "5-methylisoquinolin-3(2H)-one"),
    ("Cc1ccc2c[nH]c(=O)cc2c1",      "6-methylisoquinolin-3(2H)-one"),
    ("Cc1ccc2cc(=O)[nH]cc2c1",      "7-methylisoquinolin-3(2H)-one"),
    ("Cc1cccc2cc(=O)[nH]cc12",      "8-methylisoquinolin-3(2H)-one"),
    # cinnolin-4(1H)-one (CH at 3,5,6,7,8)
    ("O=c1cn[nH]c2ccccc12",         "cinnolin-4(1H)-one"),
    ("Cc1n[nH]c2ccccc2c1=O",        "3-methylcinnolin-4(1H)-one"),
    ("Cc1cccc2[nH]ncc(=O)c12",      "5-methylcinnolin-4(1H)-one"),
    ("Cc1ccc2[nH]ncc(=O)c2c1",      "6-methylcinnolin-4(1H)-one"),
    ("Cc1ccc2c(=O)cn[nH]c2c1",      "7-methylcinnolin-4(1H)-one"),
    ("Cc1cccc2c(=O)cn[nH]c12",      "8-methylcinnolin-4(1H)-one"),
    # phthalazin-1(2H)-one (CH at 4,5,6,7,8)
    ("O=c1[nH]ncc2ccccc12",         "phthalazin-1(2H)-one"),
    ("Cc1n[nH]c(=O)c2ccccc12",      "4-methylphthalazin-1(2H)-one"),
    ("Cc1cccc2c(=O)[nH]ncc12",      "5-methylphthalazin-1(2H)-one"),
    ("Cc1ccc2c(=O)[nH]ncc2c1",      "6-methylphthalazin-1(2H)-one"),
    ("Cc1ccc2cn[nH]c(=O)c2c1",      "7-methylphthalazin-1(2H)-one"),
    ("Cc1cccc2cn[nH]c(=O)c12",      "8-methylphthalazin-1(2H)-one"),
    # quinazolin-4(3H)-one (CH at 2,5,6,7,8)
    ("O=c1[nH]cnc2ccccc12",         "quinazolin-4(3H)-one"),
    ("Cc1nc2ccccc2c(=O)[nH]1",      "2-methylquinazolin-4(3H)-one"),
    ("Cc1cccc2nc[nH]c(=O)c12",      "5-methylquinazolin-4(3H)-one"),
    ("Cc1ccc2nc[nH]c(=O)c2c1",      "6-methylquinazolin-4(3H)-one"),
    ("Cc1ccc2c(=O)[nH]cnc2c1",      "7-methylquinazolin-4(3H)-one"),
    ("Cc1cccc2c(=O)[nH]cnc12",      "8-methylquinazolin-4(3H)-one"),
    # quinazoline-2,4(1H,3H)-dione (CH at 5,6,7,8)
    ("O=c1[nH]c(=O)c2ccccc2[nH]1", "quinazoline-2,4(1H,3H)-dione"),
    ("Cc1cccc2[nH]c(=O)[nH]c(=O)c12", "5-methylquinazoline-2,4(1H,3H)-dione"),
    ("Cc1ccc2[nH]c(=O)[nH]c(=O)c2c1", "6-methylquinazoline-2,4(1H,3H)-dione"),
    ("Cc1ccc2c(=O)[nH]c(=O)[nH]c2c1", "7-methylquinazoline-2,4(1H,3H)-dione"),
    ("Cc1cccc2c(=O)[nH]c(=O)[nH]c12", "8-methylquinazoline-2,4(1H,3H)-dione"),
    # quinoxalin-2(1H)-one (CH at 3,5,6,7,8)
    ("O=c1cnc2ccccc2[nH]1",         "quinoxalin-2(1H)-one"),
    ("Cc1nc2ccccc2[nH]c1=O",        "3-methylquinoxalin-2(1H)-one"),
    ("Cc1cccc2[nH]c(=O)cnc12",      "5-methylquinoxalin-2(1H)-one"),
    ("Cc1ccc2[nH]c(=O)cnc2c1",      "6-methylquinoxalin-2(1H)-one"),
    ("Cc1ccc2ncc(=O)[nH]c2c1",      "7-methylquinoxalin-2(1H)-one"),
    ("Cc1cccc2ncc(=O)[nH]c12",      "8-methylquinoxalin-2(1H)-one"),
])
def test_phase719(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
