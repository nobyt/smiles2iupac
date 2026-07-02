"""Phase 768: benzo[g]cinnoline and benzo[h]cinnoline γ-ol/thiol → tautomers (IUPAC 2013).

Analogous to cinnoline-4(1H)-one (gamma position, not the direct alpha):
- benzo[g]cinnoline: C5 (gamma to N7, outer N8 gets H) → benzo[g]cinnolin-5(8H)-one
- benzo[h]cinnoline: C7 (gamma to N9, outer N10 gets H) → benzo[h]cinnolin-7(10H)-one
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[g]cinnoline
    ("Oc1cnnc2cc3ccccc3cc12",          "benzo[g]cinnolin-5(8H)-one"),
    ("Sc1cnnc2cc3ccccc3cc12",          "benzo[g]cinnolin-5(8H)-thione"),
    # benzo[h]cinnoline
    ("Oc1cnnc2c1ccc1ccccc12",          "benzo[h]cinnolin-7(10H)-one"),
    ("Sc1cnnc2c1ccc1ccccc12",          "benzo[h]cinnolin-7(10H)-thione"),
    # Regression: parent rings unaffected
    ("c1ccc2cc3nnccc3cc2c1",           "benzo[g]cinnoline"),
    ("c1ccc2c(c1)ccc1ccnnc12",         "benzo[h]cinnoline"),
    # Regression: direct alpha positions not converted
    ("Oc1cc2cc3ccccc3cc2nn1",          "benzo[g]cinnolin-6-ol"),
    ("Oc1cc2ccc3ccccc3c2nn1",          "benzo[h]cinnolin-8-ol"),
    # Regression: Phase 767 unchanged
    ("Oc1cnc2ccc3ccccc3c2n1",          "benzo[f]quinoxalin-3(2H)-one"),
])
def test_phase768_benzo_cinnoline_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
