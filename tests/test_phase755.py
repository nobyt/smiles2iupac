"""Phase 755: 1,3-thiazol-2-ol → thiazolinone; 9H-purin-6-ol → purinone (IUPAC 2013).

Extends 5-membered ring tautomer rules (Phase 754) to 1,3-thiazole and purine:
- 1,3-thiazol-2-ol → 1,3-thiazol-2(3H)-one  (and 2-thiol → 2(3H)-thione)
- 9H-purin-6-ol → 9H-purin-6(1H)-one  (hypoxanthine; also 2-amino → guanine)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,3-thiazol-2-ol → 1,3-thiazol-2(3H)-one
    ("Oc1nccs1",                    "1,3-thiazol-2(3H)-one"),
    ("Sc1nccs1",                    "1,3-thiazol-2(3H)-thione"),
    # 9H-purin-6-ol → 9H-purin-6(1H)-one (hypoxanthine)
    ("Oc1ncnc2[nH]cnc12",          "9H-purin-6(1H)-one"),
    # 2-amino-9H-purin-6-ol → 2-amino-9H-purin-6(1H)-one (guanine)
    ("Oc1nc(N)nc2[nH]cnc12",       "2-amino-9H-purin-6(1H)-one"),
    # Regression: keto form unchanged
    ("O=C1NC=CS1",                  "1,3-thiazol-2(3H)-one"),
    # Regression: parent rings unaffected
    ("c1nccs1",                     "1,3-thiazole"),
    ("c1ncc2nc[nH]c2n1",           "9H-purine"),
    # Regression: Phase 754 cases unchanged
    ("Oc1ncc[nH]1",                "1H-imidazol-2(3H)-one"),
    ("Oc1ncco1",                   "1,3-oxazol-2(3H)-one"),
])
def test_phase755_thiazole_and_purine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
