"""Phase 662: 3,4-dihydroquinolin-2(1H)-one and 3,4-dihydroisoquinolin-1(2H)-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 3,4-dihydroquinolin-2(1H)-one: N(1)-C(2,=O)-C(3)-C(4)-C4a-...-C8a
    # (parent and 1-methyl in Phase 409; C2=O not methylable)
    ("O=C1CCc2ccccc2N1",    "3,4-dihydroquinolin-2(1H)-one"),
    ("O=C1C(C)Cc2ccccc2N1", "3-methyl-3,4-dihydroquinolin-2(1H)-one"),
    ("O=C1CC(C)c2ccccc2N1", "4-methyl-3,4-dihydroquinolin-2(1H)-one"),
    ("O=C1CCc2c(C)cccc2N1", "5-methyl-3,4-dihydroquinolin-2(1H)-one"),
    ("O=C1CCc2cc(C)ccc2N1", "6-methyl-3,4-dihydroquinolin-2(1H)-one"),
    ("O=C1CCc2ccc(C)cc2N1", "7-methyl-3,4-dihydroquinolin-2(1H)-one"),
    ("O=C1CCc2cccc(C)c2N1", "8-methyl-3,4-dihydroquinolin-2(1H)-one"),
    # 3,4-dihydroisoquinolin-1(2H)-one: C(1,=O)-N(2)-C(3)-C(4)-C4a-...-C8a
    # (parent and 2-methyl in Phase 409; C1=O not methylable)
    ("O=C1NCCc2ccccc21",    "3,4-dihydroisoquinolin-1(2H)-one"),
    ("O=C1NC(C)Cc2ccccc21", "3-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    ("O=C1NCC(C)c2ccccc21", "4-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    ("O=C1NCCc2c(C)cccc21", "5-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    ("O=C1NCCc2cc(C)ccc21", "6-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    ("O=C1NCCc2ccc(C)cc21", "7-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    ("O=C1NCCc2cccc(C)c21", "8-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
])
def test_phase662(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
