"""Phase 669: 1H-indazol-3(2H)-one and phenanthridin-6(5H)-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-indazol-3(2H)-one: C(3,=O)-N(2,H)-N(1,H)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 417; C3=O, N2-H, N1-H not methylable)
    ("O=c1[nH][nH]c2ccccc12",    "1H-indazol-3(2H)-one"),
    ("O=c1[nH][nH]c2c(C)cccc12", "4-methyl-1H-indazol-3(2H)-one"),
    ("O=c1[nH][nH]c2cc(C)ccc12", "5-methyl-1H-indazol-3(2H)-one"),
    ("O=c1[nH][nH]c2ccc(C)cc12", "6-methyl-1H-indazol-3(2H)-one"),
    ("O=c1[nH][nH]c2cccc(C)c12", "7-methyl-1H-indazol-3(2H)-one"),
    # phenanthridin-6(5H)-one: C(6,=O)-N(5,H)-C4b-C1..C4-C4a / C10a-C7..C10-C6a
    # (parent in Phase 417; C6=O and N5-H not methylable; 8 unique positions)
    ("O=c1[nH]c2ccccc2c2ccccc12",    "phenanthridin-6(5H)-one"),
    ("O=c1[nH]c2cccc(C)c2c2ccccc12", "1-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2ccc(C)cc2c2ccccc12", "2-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2cc(C)ccc2c2ccccc12", "3-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2c(C)cccc2c2ccccc12", "4-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2ccccc2c2cccc(C)c12", "7-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2ccccc2c2ccc(C)cc12", "8-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2ccccc2c2cc(C)ccc12", "9-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2ccccc2c2c(C)cccc12", "10-methylphenanthridin-6(5H)-one"),
])
def test_phase669(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
