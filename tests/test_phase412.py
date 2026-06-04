"""Phase 412: Benzimidazol-2-one, benzothiazol-2-one, benzoxazol-2-one retained names.

IUPAC 2013 P-31.1.3: benzo-fused 5-membered lactam/lactone C2=O derivatives
of benzimidazole, benzothiazole, and benzoxazole.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-benzimidazol-2(3H)-one (both NH groups)
    ("O=c1[nH]c2ccccc2[nH]1",         "1H-benzimidazol-2(3H)-one"),
    # 1,3-benzothiazol-2(3H)-one (S at 1, NH at 3)
    ("O=c1[nH]c2ccccc2s1",             "1,3-benzothiazol-2(3H)-one"),
    # 1,3-benzoxazol-2(3H)-one (O at 1, NH at 3)
    ("O=c1[nH]c2ccccc2o1",             "1,3-benzoxazol-2(3H)-one"),
    # regression: 1H-benzimidazole unchanged
    ("c1ccc2[nH]cnc2c1",               "1H-benzimidazole"),
    # regression: 1,3-benzothiazole unchanged
    ("c1ccc2scnc2c1",                   "1,3-benzothiazole"),
    # regression: 1,3-benzoxazole unchanged
    ("c1ccc2ocnc2c1",                   "1,3-benzoxazole"),
    # regression: phthalimide unchanged
    ("O=C1NC(=O)c2ccccc21",            "phthalimide"),
    # regression: benzene unchanged
    ("c1ccccc1",                        "benzene"),
])
def test_phase412_benzazolones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
