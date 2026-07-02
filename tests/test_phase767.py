"""Phase 767: benzo[f]quinoxaline, benzo[g/f]phthalazine α-ol/thiol → tautomers (IUPAC 2013).

- benzo[f]quinoxaline: positions 3(N2), 4(N5)
- benzo[g]phthalazine: position 5(N6) [symmetric ring]
- benzo[f]phthalazine: positions 5(N4), 2(N3)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[f]quinoxaline
    ("Oc1cnc2ccc3ccccc3c2n1",          "benzo[f]quinoxalin-3(2H)-one"),
    ("Sc1cnc2ccc3ccccc3c2n1",          "benzo[f]quinoxalin-3(2H)-thione"),
    ("Oc1cnc2c(ccc3ccccc32)n1",        "benzo[f]quinoxalin-4(5H)-one"),
    ("Sc1cnc2c(ccc3ccccc32)n1",        "benzo[f]quinoxalin-4(5H)-thione"),
    # benzo[g]phthalazine (symmetric: C5 and C8 give the same canonical SMILES)
    ("Oc1nncc2cc3ccccc3cc12",          "benzo[g]phthalazin-5(6H)-one"),
    ("Sc1nncc2cc3ccccc3cc12",          "benzo[g]phthalazin-5(6H)-thione"),
    # benzo[f]phthalazine
    ("Oc1nncc2c1ccc1ccccc12",          "benzo[f]phthalazin-5(4H)-one"),
    ("Sc1nncc2c1ccc1ccccc12",          "benzo[f]phthalazin-5(4H)-thione"),
    ("Oc1nncc2ccc3ccccc3c12",          "benzo[f]phthalazin-2(3H)-one"),
    ("Sc1nncc2ccc3ccccc3c12",          "benzo[f]phthalazin-2(3H)-thione"),
    # Regression: parent rings unaffected
    ("c1ccc2c(c1)ccc1nccnc12",         "benzo[f]quinoxaline"),
    ("c1ccc2cc3cnncc3cc2c1",           "benzo[g]phthalazine"),
    ("c1ccc2c(c1)ccc1cnncc12",         "benzo[f]phthalazine"),
    # Regression: Phase 766 unchanged
    ("Oc1cnc2cc3ccccc3cc2n1",          "benzo[g]quinoxalin-6(5H)-one"),
])
def test_phase767_benzo_quinoxaline_phthalazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
