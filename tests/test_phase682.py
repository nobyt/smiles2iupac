"""Phase 682: phthalazin-1(2H)-one, cinnolin-4(1H)-one, quinazolin-4(3H)-one, quinolin-2(1H)-one, and quinazoline-2,4(1H,3H)-dione methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phthalazin-1(2H)-one: C(1,=O)-N(2,H)-N(3)-C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent in Phase 414; C1=O and N3 not methylable; C4 and C5-C8 methylable)
    ("O=c1[nH]ncc2ccccc12",    "phthalazin-1(2H)-one"),
    ("O=c1[nH]nc(C)c2ccccc12", "4-methylphthalazin-1(2H)-one"),
    ("O=c1[nH]ncc2c(C)cccc12", "5-methylphthalazin-1(2H)-one"),
    ("O=c1[nH]ncc2cc(C)ccc12", "6-methylphthalazin-1(2H)-one"),
    ("O=c1[nH]ncc2ccc(C)cc12", "7-methylphthalazin-1(2H)-one"),
    ("O=c1[nH]ncc2cccc(C)c12", "8-methylphthalazin-1(2H)-one"),
    # cinnolin-4(1H)-one: C(4,=O)-C(3)-N(2)-N(1,H)-C8a-C8-C7-C6-C5-C4a
    # (parent in Phase 414; C4=O and N2 not methylable; C3 and C5-C8 methylable)
    ("O=c1cn[nH]c2ccccc12",    "cinnolin-4(1H)-one"),
    ("O=c1c(C)n[nH]c2ccccc12", "3-methylcinnolin-4(1H)-one"),
    ("O=c1cn[nH]c2c(C)cccc12", "8-methylcinnolin-4(1H)-one"),
    ("O=c1cn[nH]c2cc(C)ccc12", "7-methylcinnolin-4(1H)-one"),
    ("O=c1cn[nH]c2ccc(C)cc12", "6-methylcinnolin-4(1H)-one"),
    ("O=c1cn[nH]c2cccc(C)c12", "5-methylcinnolin-4(1H)-one"),
    # quinazolin-4(3H)-one: C(4,=O)-N(3,H)-C(2)-N(1)-C8a-C8-C7-C6-C5-C4a
    # (parent in Phase 414; C4=O and N1 not methylable; C2 and C5-C8 methylable)
    ("O=c1[nH]cnc2ccccc12",    "quinazolin-4(3H)-one"),
    ("O=c1[nH]c(C)nc2ccccc12", "2-methylquinazolin-4(3H)-one"),
    ("O=c1[nH]cnc2c(C)cccc12", "8-methylquinazolin-4(3H)-one"),
    ("O=c1[nH]cnc2cc(C)ccc12", "7-methylquinazolin-4(3H)-one"),
    ("O=c1[nH]cnc2ccc(C)cc12", "6-methylquinazolin-4(3H)-one"),
    ("O=c1[nH]cnc2cccc(C)c12", "5-methylquinazolin-4(3H)-one"),
    # quinolin-2(1H)-one: N(1,H)-C(2,=O)-C(3)-C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent in Phase 415; C2=O not methylable; C3, C4, C5-C8 methylable)
    ("O=c1ccc2ccccc2[nH]1",    "quinolin-2(1H)-one"),
    ("O=c1c(C)cc2ccccc2[nH]1", "3-methylquinolin-2(1H)-one"),
    ("O=c1cc(C)c2ccccc2[nH]1", "4-methylquinolin-2(1H)-one"),
    ("O=c1ccc2c(C)cccc2[nH]1", "5-methylquinolin-2(1H)-one"),
    ("O=c1ccc2cc(C)ccc2[nH]1", "6-methylquinolin-2(1H)-one"),
    ("O=c1ccc2ccc(C)cc2[nH]1", "7-methylquinolin-2(1H)-one"),
    ("O=c1ccc2cccc(C)c2[nH]1", "8-methylquinolin-2(1H)-one"),
    # quinazoline-2,4(1H,3H)-dione: N(1,H)-C(2,=O)-N(3,H)-C(4,=O)-C4a-C5-C6-C7-C8-C8a
    # (parent in Phase 415; both C=O not methylable; C5-C8 aromatic methylable)
    ("O=c1[nH]c(=O)c2ccccc2[nH]1",    "quinazoline-2,4(1H,3H)-dione"),
    ("O=c1[nH]c(=O)c2c(C)cccc2[nH]1", "5-methylquinazoline-2,4(1H,3H)-dione"),
    ("O=c1[nH]c(=O)c2cc(C)ccc2[nH]1", "6-methylquinazoline-2,4(1H,3H)-dione"),
    ("O=c1[nH]c(=O)c2ccc(C)cc2[nH]1", "7-methylquinazoline-2,4(1H,3H)-dione"),
    ("O=c1[nH]c(=O)c2cccc(C)c2[nH]1", "8-methylquinazoline-2,4(1H,3H)-dione"),
])
def test_phase682(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
