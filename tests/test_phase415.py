"""Phase 415: Anthracen-9(10H)-one, quinolin-2(1H)-one, quinazoline-2,4(1H,3H)-dione.

IUPAC 2013: tricyclic ketone (anthrone), bicyclic lactam, and bicyclic dione
retained/systematic names.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # anthracen-9(10H)-one (anthrone) — tricyclic ketone with CH2 bridge
    ("O=C1c2ccccc2Cc2ccccc21",          "anthracen-9(10H)-one"),
    # quinolin-2(1H)-one — benzo-fused lactam at C2
    ("O=c1ccc2ccccc2[nH]1",             "quinolin-2(1H)-one"),
    # quinazoline-2,4(1H,3H)-dione — two exo C=O groups
    ("O=c1[nH]c(=O)c2ccccc2[nH]1",     "quinazoline-2,4(1H,3H)-dione"),
    # regression: xanthen-9-one unchanged
    ("O=C1c2ccccc2Oc2ccccc21",           "xanthen-9-one"),
    # regression: anthracene unchanged
    ("c1ccc2cc3ccccc3cc2c1",             "anthracene"),
    # regression: quinoline unchanged
    ("c1ccc2ncccc2c1",                   "quinoline"),
    # regression: quinazoline unchanged
    ("c1ccc2ncncc2c1",                   "quinazoline"),
    # regression: benzene unchanged
    ("c1ccccc1",                          "benzene"),
])
def test_phase415_anthrone_quinolinone_dione(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
