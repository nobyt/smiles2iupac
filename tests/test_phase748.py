"""Phase 748: hydroxy benzo-fused diazines → preferred lactam tautomers (IUPAC 2013).

Extends Phase 744–747 to quinoxaline, quinazoline, cinnoline, and phthalazine:
  2-hydroxyquinoxaline   → quinoxalin-2(1H)-one
  2/4-hydroxyquinazoline → quinazolin-2(1H)-one / quinazolin-4(3H)-one
  4-hydroxycinnoline     → cinnolin-4(1H)-one
  1-hydroxyphthalazine   → phthalazin-1(2H)-one
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # quinoxaline: α-OH → lactam
    ("Oc1cnc2ccccc2n1",            "quinoxalin-2(1H)-one"),
    # quinazoline: 2-OH and 4-OH → lactam
    ("Oc1nc2ccccc2cn1",            "quinazolin-2(1H)-one"),
    ("Oc1ncnc2ccccc12",            "quinazolin-4(3H)-one"),
    # cinnoline: 4-OH → lactam
    ("Oc1cnnc2ccccc12",            "cinnolin-4(1H)-one"),
    # phthalazine: 1-OH → lactam
    ("Oc1nncc2ccccc12",            "phthalazin-1(2H)-one"),
    # Regression: keto SMILES unchanged
    ("O=c1cnc2ccccc2[nH]1",        "quinoxalin-2(1H)-one"),
    ("O=c1[nH]cnc2ccccc12",        "quinazolin-4(3H)-one"),
    ("O=c1cn[nH]c2ccccc12",        "cinnolin-4(1H)-one"),
    ("O=c1[nH]ncc2ccccc12",        "phthalazin-1(2H)-one"),
    # Regression: parent rings unaffected
    ("c1ccc2nccnc2c1",             "quinoxaline"),
    ("c1ccc2ncncc2c1",             "quinazoline"),
    ("c1ccc2nnccc2c1",             "cinnoline"),
    ("c1ccc2cnncc2c1",             "phthalazine"),
])
def test_phase748_benzo_diazine_ol_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
