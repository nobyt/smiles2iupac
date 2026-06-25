"""Phase 665: phthalazin-1(2H)-one, cinnolin-4(1H)-one, and quinazolin-4(3H)-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phthalazin-1(2H)-one: C(1,=O)-N(2,H)-N(3)-C(4)-C4a-...-C8a
    # (parent in Phase 414; C1=O, N2-H, and aromatic N3 not methylable)
    ("O=c1[nH]ncc2ccccc12",    "phthalazin-1(2H)-one"),
    ("O=c1[nH]nc(C)c2ccccc12", "4-methylphthalazin-1(2H)-one"),
    ("O=c1[nH]ncc2c(C)cccc12", "5-methylphthalazin-1(2H)-one"),
    ("O=c1[nH]ncc2cc(C)ccc12", "6-methylphthalazin-1(2H)-one"),
    ("O=c1[nH]ncc2ccc(C)cc12", "7-methylphthalazin-1(2H)-one"),
    ("O=c1[nH]ncc2cccc(C)c12", "8-methylphthalazin-1(2H)-one"),
    # cinnolin-4(1H)-one: C(4,=O)-C(3)-N(2)-N(1,H)-C8a-...-C4a
    # (parent in Phase 414; C4=O, N1-H, and aromatic N2 not methylable)
    ("O=c1cn[nH]c2ccccc12",    "cinnolin-4(1H)-one"),
    ("O=c1c(C)n[nH]c2ccccc12", "3-methylcinnolin-4(1H)-one"),
    ("O=c1cn[nH]c2cccc(C)c12", "5-methylcinnolin-4(1H)-one"),
    ("O=c1cn[nH]c2ccc(C)cc12", "6-methylcinnolin-4(1H)-one"),
    ("O=c1cn[nH]c2cc(C)ccc12", "7-methylcinnolin-4(1H)-one"),
    ("O=c1cn[nH]c2c(C)cccc12", "8-methylcinnolin-4(1H)-one"),
    # quinazolin-4(3H)-one: C(4,=O)-N(3,H)-C(2)-N(1)-C8a-...-C4a
    # (parent in Phase 414; C4=O, N3-H, and aromatic N1 not methylable)
    ("O=c1[nH]cnc2ccccc12",    "quinazolin-4(3H)-one"),
    ("O=c1[nH]c(C)nc2ccccc12", "2-methylquinazolin-4(3H)-one"),
    ("O=c1[nH]cnc2cccc(C)c12", "5-methylquinazolin-4(3H)-one"),
    ("O=c1[nH]cnc2ccc(C)cc12", "6-methylquinazolin-4(3H)-one"),
    ("O=c1[nH]cnc2cc(C)ccc12", "7-methylquinazolin-4(3H)-one"),
    ("O=c1[nH]cnc2c(C)cccc12", "8-methylquinazolin-4(3H)-one"),
])
def test_phase665(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
