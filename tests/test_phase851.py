"""Phase 851: N-methyl benzimidazolone/benzoxazolone/benzothiazolone → -2-one suffix.

Two fixes:
1. Locant map: N-3 (atom 6) in 1,3-benzoxazole and 1,3-benzothiazole rings
   was mapped to None, causing N-3 substituents to be silently dropped.
   Now mapped to locant 3.
2. oxo→-one conversion in _apply_hetero_suffixes for these three ring systems.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-benzimidazol-2(3H)-one → N-3 methyl drops (3H)
    ("Cn1c(=O)[nH]c2ccccc21",       "3-methyl-1H-benzimidazol-2-one"),
    # 1,3-benzoxazol-2(3H)-one → N-3 methyl drops (3H)
    ("Cn1c(=O)oc2ccccc21",          "3-methyl-1,3-benzoxazol-2-one"),
    # 1,3-benzothiazol-2(3H)-one → N-3 methyl drops (3H)
    ("Cn1c(=O)sc2ccccc21",          "3-methyl-1,3-benzothiazol-2-one"),
    # NH parents still correct (via existing hydroxy path)
    ("O=c1sc2ccccc2[nH]1",          "1,3-benzothiazol-2(3H)-one"),
    ("O=c1oc2ccccc2[nH]1",          "1,3-benzoxazol-2(3H)-one"),
    # Parent heterocycles unaffected
    ("c1ccc2ocnc2c1",               "1,3-benzoxazole"),
    ("c1ccc2scnc2c1",               "1,3-benzothiazole"),
    # C-substituents on benzoxazole/benzothiazole unaffected
    ("Cc1nc2ccccc2o1",              "2-methyl-1,3-benzoxazole"),
    ("Cc1ccc2ocnc2c1",              "5-methyl-1,3-benzoxazole"),
])
def test_phase851(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
