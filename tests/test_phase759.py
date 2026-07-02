"""Phase 759: phthalazine-1,4-diol and 1,2,4,5-tetrazin-3-ol → tautomers (IUPAC 2013).

- phthalazine-1,4-diol → phthalazine-1,4(2H,3H)-dione (and dithiol → dithione)
- 1,2,4,5-tetrazin-3-ol → 1,2,4,5-tetrazin-3(2H)-one (and 3-thiol → 3(2H)-thione)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phthalazine-1,4-diol → phthalazine-1,4(2H,3H)-dione
    ("Oc1nnc(O)c2ccccc12",         "phthalazine-1,4(2H,3H)-dione"),
    ("Sc1nnc(S)c2ccccc12",         "phthalazine-1,4(2H,3H)-dithione"),
    # 1,2,4,5-tetrazin-3-ol → 1,2,4,5-tetrazin-3(2H)-one
    ("Oc1nncnn1",                  "1,2,4,5-tetrazin-3(2H)-one"),
    ("Sc1nncnn1",                  "1,2,4,5-tetrazin-3(2H)-thione"),
    # Regression: mono-OH phthalazine unchanged (Phase 748)
    ("Oc1cnnc2ccccc21",            "cinnolin-4(1H)-one"),
    # Regression: parent rings unaffected
    ("c1ccc2cnncc2c1",             "phthalazine"),
    ("c1nncnn1",                   "1,2,4,5-tetrazine"),
    # Regression: Phase 758 unchanged
    ("Oc1nc(O)c2ccccc2n1",         "quinazoline-2,4(1H,3H)-dione"),
    ("Oc1nc2ccccc2nc1O",           "quinoxaline-2,3(1H,4H)-dione"),
])
def test_phase759_phthalazine_tetrazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
