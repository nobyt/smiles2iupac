"""Phase 769: benzo[f/h/g]quinoline additional α-ol/thiol → tautomers (IUPAC 2013).

- benzo[f]quinoline: C4 (alpha to N5) → benzo[f]quinolin-4(5H)-one
- benzo[h]quinoline: C9 (alpha to N10) → benzo[h]quinolin-9(10H)-one
- benzo[g]quinoline: C2 (alpha to N1) → benzo[g]quinolin-2(1H)-one
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[f]quinoline C4
    ("Oc1ccc2c(ccc3ccccc32)n1",         "benzo[f]quinolin-4(5H)-one"),
    ("Sc1ccc2c(ccc3ccccc32)n1",         "benzo[f]quinolin-4(5H)-thione"),
    # benzo[h]quinoline C9
    ("Oc1ccc2ccc3ccccc3c2n1",           "benzo[h]quinolin-9(10H)-one"),
    ("Sc1ccc2ccc3ccccc3c2n1",           "benzo[h]quinolin-9(10H)-thione"),
    # benzo[g]quinoline C2
    ("Oc1ccc2cc3ccccc3cc2n1",           "benzo[g]quinolin-2(1H)-one"),
    ("Sc1ccc2cc3ccccc3cc2n1",           "benzo[g]quinolin-2(1H)-thione"),
    # Regression: parent rings unaffected
    ("c1ccc2c(c1)ccc1ncccc12",          "benzo[f]quinoline"),
    ("c1ccc2c(c1)ccc1cccnc12",          "benzo[h]quinoline"),
    ("c1ccc2cc3ccccc3cc2n1",            "benzo[g]quinoline"),
    # Regression: Phase 761/762 unchanged
    ("Oc1ccnc2ccc3ccccc3c12",           "benzo[f]quinolin-2(1H)-one"),
    ("Oc1cccc2ccc3cccnc3c21",           "benzo[h]quinolin-1(10H)-one"),
])
def test_phase769_benzo_quinoline_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
