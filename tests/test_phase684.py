"""Phase 684: 1H-indazol-3(2H)-one, phenanthridin-6(5H)-one, and phenanthridine methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-indazol-3(2H)-one: N(1,H)-N(2,H)-C(3,=O)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 417; C3=O not methylable; C4-C7 aromatic methylable)
    ("O=c1[nH][nH]c2ccccc12",    "1H-indazol-3(2H)-one"),
    ("O=c1[nH][nH]c2c(C)cccc12", "4-methyl-1H-indazol-3(2H)-one"),
    ("O=c1[nH][nH]c2cc(C)ccc12", "5-methyl-1H-indazol-3(2H)-one"),
    ("O=c1[nH][nH]c2ccc(C)cc12", "6-methyl-1H-indazol-3(2H)-one"),
    ("O=c1[nH][nH]c2cccc(C)c12", "7-methyl-1H-indazol-3(2H)-one"),
    # phenanthridin-6(5H)-one: C(6,=O)-N(5,H)-C4b-C4-C3-C2-C1-C4a-C10a-C10-C9-C8-C7-C6a
    # (parent in Phase 417; C6=O not methylable; C1-C4 and C7-C10 aromatic methylable)
    ("O=c1[nH]c2ccccc2c2ccccc12",    "phenanthridin-6(5H)-one"),
    ("O=c1[nH]c2cccc(C)c2c2ccccc12", "1-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2ccc(C)cc2c2ccccc12", "2-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2cc(C)ccc2c2ccccc12", "3-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2c(C)cccc2c2ccccc12", "4-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2ccccc2c2cccc(C)c12", "7-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2ccccc2c2ccc(C)cc12", "8-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2ccccc2c2cc(C)ccc12", "9-methylphenanthridin-6(5H)-one"),
    ("O=c1[nH]c2ccccc2c2c(C)cccc12", "10-methylphenanthridin-6(5H)-one"),
    # phenanthridine: N at 5; C1-C4, C6-C10 all aromatic methylable
    ("c1ccc2c(c1)cnc1ccccc12",    "phenanthridine"),
    ("Cc1cccc2cnc3ccccc3c12",     "1-methylphenanthridine"),
    ("Cc1ccc2cnc3ccccc3c2c1",     "2-methylphenanthridine"),
    ("Cc1ccc2c(c1)cnc1ccccc12",   "3-methylphenanthridine"),
    ("Cc1cccc2c1cnc1ccccc12",     "4-methylphenanthridine"),
    ("c1ccc2c(c1)c(C)nc1ccccc12", "6-methylphenanthridine"),
    ("c1ccc2c(c1)cnc1c(C)cccc12", "7-methylphenanthridine"),
    ("c1ccc2c(c1)cnc1cc(C)ccc12", "8-methylphenanthridine"),
    ("c1ccc2c(c1)cnc1ccc(C)cc12", "9-methylphenanthridine"),
    ("c1ccc2c(c1)cnc1cccc(C)c12", "10-methylphenanthridine"),
])
def test_phase684(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
