"""Phase 663: quinolin-2(1H)-one and quinolin-4(1H)-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # quinolin-2(1H)-one: N(1,H)-C(2,=O)-C(3)-C(4)-C4a-...-C8a
    # (parent in Phase 415; C2=O and N1-H not methylable)
    ("O=c1ccc2ccccc2[nH]1",    "quinolin-2(1H)-one"),
    ("O=c1c(C)cc2ccccc2[nH]1", "3-methylquinolin-2(1H)-one"),
    ("O=c1cc(C)c2ccccc2[nH]1", "4-methylquinolin-2(1H)-one"),
    ("O=c1ccc2c(C)cccc2[nH]1", "5-methylquinolin-2(1H)-one"),
    ("O=c1ccc2cc(C)ccc2[nH]1", "6-methylquinolin-2(1H)-one"),
    ("O=c1ccc2ccc(C)cc2[nH]1", "7-methylquinolin-2(1H)-one"),
    ("O=c1ccc2cccc(C)c2[nH]1", "8-methylquinolin-2(1H)-one"),
    # quinolin-4(1H)-one: N(1,H)-C(2)-C(3)-C(4,=O)-C4a-...-C8a
    # (parent in Phase 416; C4=O and N1-H not methylable)
    ("O=c1cc[nH]c2ccccc12",    "quinolin-4(1H)-one"),
    ("O=c1cc(C)[nH]c2ccccc12", "2-methylquinolin-4(1H)-one"),
    ("O=c1c(C)c[nH]c2ccccc12", "3-methylquinolin-4(1H)-one"),
    ("O=c1cc[nH]c2cccc(C)c12", "5-methylquinolin-4(1H)-one"),
    ("O=c1cc[nH]c2ccc(C)cc12", "6-methylquinolin-4(1H)-one"),
    ("O=c1cc[nH]c2cc(C)ccc12", "7-methylquinolin-4(1H)-one"),
    ("O=c1cc[nH]c2c(C)cccc12", "8-methylquinolin-4(1H)-one"),
])
def test_phase663(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
