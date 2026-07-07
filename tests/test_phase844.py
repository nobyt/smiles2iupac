"""Phase 844: drop (nH) inline indicated-H when the N at that locant is substituted.

IUPAC 2013 P-14.7: indicated H notation (nH) is used only for H atoms at sp3
positions. When the N instead carries a substituent, the (nH) must be omitted.

This extends Phase 843 to retained names with embedded (nH) patterns, via a
general regex substitution that strips (nH) from the base name whenever locant
n appears in the substituent list.

Affected patterns:
- 3,4-dihydroquinolin-2(1H)-one: drop (1H) when N1 is substituted.
- 3,4-dihydroisoquinolin-1(2H)-one: drop (2H) when N2 is substituted.
- 2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one: drop (4H) when N4 is substituted.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 3,4-dihydroquinolin-2(1H)-one: N1-H → keep (1H)
    ("O=C1CCc2ccccc2N1",                 "3,4-dihydroquinolin-2(1H)-one"),
    # N1-methyl → drop (1H)
    ("CN1C(=O)CCc2ccccc21",              "1-methyl-3,4-dihydroquinolin-2-one"),
    # C3-methyl: N1 still has H → keep (1H)
    ("O=C1C(C)Cc2ccccc2N1",             "3-methyl-3,4-dihydroquinolin-2(1H)-one"),
    # 3,4-dihydroisoquinolin-1(2H)-one: N2-H → keep (2H)
    ("O=C1NCCc2ccccc21",                 "3,4-dihydroisoquinolin-1(2H)-one"),
    # N2-methyl → drop (2H)
    ("CN1CCc2ccccc2C1=O",               "2-methyl-3,4-dihydroisoquinolin-1-one"),
    # C3-methyl: N2 still has H → keep (2H)
    ("O=C1NC(C)Cc2ccccc21",             "3-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    # 2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one: N4-H → keep (4H)
    ("O=C1NCCNc2ccccc21",               "2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one"),
    # N4-methyl → drop (4H); N1 still has H so 1H- stays
    ("CN1CCNc2ccccc2C1=O",              "4-methyl-2,3-dihydro-1H-1,4-benzodiazepin-5-one"),
])
def test_phase844(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
