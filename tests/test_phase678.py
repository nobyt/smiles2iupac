"""Phase 678: 3,4-dihydroquinolin-2(1H)-one and 3,4-dihydroisoquinolin-1(2H)-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 3,4-dihydroquinolin-2(1H)-one: N(1,H)-C(2,=O)-C(3)-C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent in Phase 409; C2=O not methylable; N1-H and C3-C8 methylable)
    ("O=C1CCc2ccccc2N1",    "3,4-dihydroquinolin-2(1H)-one"),
    ("CN1C(=O)CCc2ccccc21", "1-methyl-3,4-dihydroquinolin-2(1H)-one"),
    ("O=C1C(C)Cc2ccccc2N1", "3-methyl-3,4-dihydroquinolin-2(1H)-one"),
    ("O=C1CC(C)c2ccccc2N1", "4-methyl-3,4-dihydroquinolin-2(1H)-one"),
    ("O=C1CCc2c(C)cccc2N1", "5-methyl-3,4-dihydroquinolin-2(1H)-one"),
    ("O=C1CCc2cc(C)ccc2N1", "6-methyl-3,4-dihydroquinolin-2(1H)-one"),
    ("O=C1CCc2ccc(C)cc2N1", "7-methyl-3,4-dihydroquinolin-2(1H)-one"),
    ("O=C1CCc2cccc(C)c2N1", "8-methyl-3,4-dihydroquinolin-2(1H)-one"),
    # 3,4-dihydroisoquinolin-1(2H)-one: C(1,=O)-N(2,H)-C(3)-C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent in Phase 409; C1=O not methylable; N2-H and C3-C8 methylable)
    ("O=C1NCCc2ccccc21",    "3,4-dihydroisoquinolin-1(2H)-one"),
    ("CN1CCc2ccccc2C1=O",   "2-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    ("O=C1NC(C)Cc2ccccc21", "3-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    ("O=C1NCC(C)c2ccccc21", "4-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    ("O=C1NCCc2c(C)cccc21", "5-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    ("O=C1NCCc2cc(C)ccc21", "6-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    ("O=C1NCCc2ccc(C)cc21", "7-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    ("O=C1NCCc2cccc(C)c21", "8-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
])
def test_phase678(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
