"""Phase 762: benzo[h]quinoline-1-ol → benzo[h]quinolin-1(10H)-one (IUPAC 2013).

Extends tricyclic tautomeric conversions (Phase 761) to benzo[h]quinoline:
- benzo[h]quinoline-1-ol → benzo[h]quinolin-1(10H)-one
- benzo[h]quinoline-1-thiol → benzo[h]quinolin-1(10H)-thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1cccc2ccc3cccnc3c21",          "benzo[h]quinolin-1(10H)-one"),
    ("Sc1cccc2ccc3cccnc3c21",          "benzo[h]quinolin-1(10H)-thione"),
    # Regression: parent ring unaffected
    ("c1ccc2c(c1)ccc1cccnc12",         "benzo[h]quinoline"),
    # Regression: Phase 761 unchanged
    ("Oc1ccnc2ccc3ccccc3c12",          "benzo[f]quinolin-2(1H)-one"),
    ("Oc1c2ccccc2nc2ccccc12",          "acridin-9(10H)-one"),
])
def test_phase762_benzohquinoline_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
