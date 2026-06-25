"""Phase 683: quinolin-4(1H)-one, isoquinolin-1(2H)-one, and isoquinolin-3(2H)-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # quinolin-4(1H)-one: C(4,=O)-C(3)-C(2)-N(1,H)-C8a-C8-C7-C6-C5-C4a
    # (parent in Phase 416; C4=O not methylable; C2, C3, C5-C8 methylable)
    ("O=c1cc[nH]c2ccccc12",    "quinolin-4(1H)-one"),
    ("O=c1cc(C)[nH]c2ccccc12", "2-methylquinolin-4(1H)-one"),
    ("O=c1c(C)c[nH]c2ccccc12", "3-methylquinolin-4(1H)-one"),
    ("O=c1cc[nH]c2c(C)cccc12", "8-methylquinolin-4(1H)-one"),
    ("O=c1cc[nH]c2cc(C)ccc12", "7-methylquinolin-4(1H)-one"),
    ("O=c1cc[nH]c2ccc(C)cc12", "6-methylquinolin-4(1H)-one"),
    ("O=c1cc[nH]c2cccc(C)c12", "5-methylquinolin-4(1H)-one"),
    # isoquinolin-1(2H)-one: C(1,=O)-N(2,H)-C(3)-C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent in Phase 416; C1=O not methylable; C3, C4, C5-C8 methylable)
    ("O=c1[nH]ccc2ccccc12",    "isoquinolin-1(2H)-one"),
    ("O=c1[nH]c(C)cc2ccccc12", "3-methylisoquinolin-1(2H)-one"),
    ("O=c1[nH]cc(C)c2ccccc12", "4-methylisoquinolin-1(2H)-one"),
    ("O=c1[nH]ccc2c(C)cccc12", "5-methylisoquinolin-1(2H)-one"),
    ("O=c1[nH]ccc2cc(C)ccc12", "6-methylisoquinolin-1(2H)-one"),
    ("O=c1[nH]ccc2ccc(C)cc12", "7-methylisoquinolin-1(2H)-one"),
    ("O=c1[nH]ccc2cccc(C)c12", "8-methylisoquinolin-1(2H)-one"),
    # isoquinolin-3(2H)-one: C(3,=O)-C(4)-C4a-C5-C6-C7-C8-C8a-C(1)-N(2,H)
    # (parent in Phase 416; C3=O not methylable; C1, C4, C5-C8 methylable)
    ("O=c1cc2ccccc2c[nH]1",    "isoquinolin-3(2H)-one"),
    ("O=c1cc2ccccc2c(C)[nH]1", "1-methylisoquinolin-3(2H)-one"),
    ("O=c1c(C)c2ccccc2c[nH]1", "4-methylisoquinolin-3(2H)-one"),
    ("O=c1cc2c(C)cccc2c[nH]1", "5-methylisoquinolin-3(2H)-one"),
    ("O=c1cc2cc(C)ccc2c[nH]1", "6-methylisoquinolin-3(2H)-one"),
    ("O=c1cc2ccc(C)cc2c[nH]1", "7-methylisoquinolin-3(2H)-one"),
    ("O=c1cc2cccc(C)c2c[nH]1", "8-methylisoquinolin-3(2H)-one"),
])
def test_phase683(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
