"""Phase 745: 2/4-hydroxyquinoline and 1/3-hydroxyisoquinoline → lactam tautomers
(IUPAC 2013 preferred tautomers, extending Phase 744 to benzo-fused systems).

α/γ-hydroxyquinolines prefer quinolin-N(1H)-one; α-hydroxyisoquinolines
prefer isoquinolin-N(2H)-one. β-hydroxy positions are unaffected.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # quinoline: α (2) and γ (4) positions → lactam
    ("Oc1nc2ccccc2cc1",           "quinolin-2(1H)-one"),
    ("Oc1ccnc2ccccc12",           "quinolin-4(1H)-one"),
    # isoquinoline: α positions (1, 3 — both adjacent to N2) → lactam
    ("Oc1nccc2ccccc12",           "isoquinolin-1(2H)-one"),
    ("Oc1cc2ccccc2cn1",           "isoquinolin-3(2H)-one"),
    # Substituted: substituent prefix before the lactam base
    ("Cc1ccc2nc(O)ccc2c1",        "6-methylquinolin-2(1H)-one"),
    ("Cc1ccc2ccnc(O)c2c1",        "7-methylisoquinolin-1(2H)-one"),
    # Regression: β-hydroxy positions unaffected
    ("Oc1ccc2ncccc2c1",           "quinolin-6-ol"),
    ("Oc1cccc2ncccc12",           "quinolin-5-ol"),
    # Regression: pyridine/pyrimidine tautomers unchanged
    ("Oc1ccccn1",                 "1H-pyridin-2-one"),
    ("Oc1ccncc1",                 "1H-pyridin-4-one"),
    # Regression: quinoline keto SMILES unchanged
    ("O=c1ccc2ccccc2[nH]1",       "quinolin-2(1H)-one"),
    ("O=c1cc[nH]c2ccccc12",       "quinolin-4(1H)-one"),
])
def test_phase745_hydroxyquinoline_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
