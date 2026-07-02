"""Phase 761: acridine-9-ol, phenanthridine-6-ol, benzo[f]quinoline-2-ol → tautomers (IUPAC 2013).

Extends tautomeric conversions to tricyclic N-heterocycles:
- acridine-9-ol → acridin-9(10H)-one  (acridone)
- phenanthridine-6-ol → phenanthridin-6(5H)-one
- benzo[f]quinoline-2-ol → benzo[f]quinolin-2(1H)-one
(and thiol → thione counterparts)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acridine-9-ol → acridin-9(10H)-one
    ("Oc1c2ccccc2nc2ccccc12",          "acridin-9(10H)-one"),
    ("Sc1c2ccccc2nc2ccccc12",          "acridin-9(10H)-thione"),
    # phenanthridine-6-ol → phenanthridin-6(5H)-one
    ("Oc1nc2ccccc2c2ccccc12",          "phenanthridin-6(5H)-one"),
    ("Sc1nc2ccccc2c2ccccc12",          "phenanthridin-6(5H)-thione"),
    # benzo[f]quinoline-2-ol → benzo[f]quinolin-2(1H)-one
    ("Oc1ccnc2ccc3ccccc3c12",          "benzo[f]quinolin-2(1H)-one"),
    ("Sc1ccnc2ccc3ccccc3c12",          "benzo[f]quinolin-2(1H)-thione"),
    # Regression: parent rings unaffected
    ("c1ccc2nc3ccccc3cc2c1",           "acridine"),
    ("c1ccc2c(c1)cnc1ccccc12",         "phenanthridine"),
    ("c1ccc2c(c1)ccc1ncccc12",         "benzo[f]quinoline"),
    # Regression: Phase 760 benzotriazine unchanged
    ("Oc1nnc2ccccc2n1",                "1,2,4-benzotriazin-3(2H)-one"),
])
def test_phase761_tricyclic_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
