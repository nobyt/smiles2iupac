"""Phase 664: isoquinolin-1(2H)-one and isoquinolin-3(2H)-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isoquinolin-1(2H)-one: C(1,=O)-N(2,H)-C(3)-C(4)-C4a-...-C8a
    # (parent in Phase 416; C1=O and N2-H not methylable)
    ("O=c1[nH]ccc2ccccc12",    "isoquinolin-1(2H)-one"),
    ("O=c1[nH]c(C)cc2ccccc12", "3-methylisoquinolin-1(2H)-one"),
    ("O=c1[nH]cc(C)c2ccccc12", "4-methylisoquinolin-1(2H)-one"),
    ("O=c1[nH]ccc2c(C)cccc12", "5-methylisoquinolin-1(2H)-one"),
    ("O=c1[nH]ccc2cc(C)ccc12", "6-methylisoquinolin-1(2H)-one"),
    ("O=c1[nH]ccc2ccc(C)cc12", "7-methylisoquinolin-1(2H)-one"),
    ("O=c1[nH]ccc2cccc(C)c12", "8-methylisoquinolin-1(2H)-one"),
    # isoquinolin-3(2H)-one: C(3,=O)-C(4)-C4a-...-C8a-C(1)-N(2,H)
    # (parent in Phase 416; C3=O and N2-H not methylable)
    ("O=c1cc2ccccc2c[nH]1",    "isoquinolin-3(2H)-one"),
    ("O=c1cc2ccccc2c(C)[nH]1", "1-methylisoquinolin-3(2H)-one"),
    ("O=c1c(C)c2ccccc2c[nH]1", "4-methylisoquinolin-3(2H)-one"),
    ("O=c1cc2c(C)cccc2c[nH]1", "5-methylisoquinolin-3(2H)-one"),
    ("O=c1cc2cc(C)ccc2c[nH]1", "6-methylisoquinolin-3(2H)-one"),
    ("O=c1cc2ccc(C)cc2c[nH]1", "7-methylisoquinolin-3(2H)-one"),
    ("O=c1cc2cccc(C)c2c[nH]1", "8-methylisoquinolin-3(2H)-one"),
])
def test_phase664(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
