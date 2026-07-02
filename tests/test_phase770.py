"""Phase 770: benzo[f]cinnoline γ-ol/thiol → tautomers (IUPAC 2013).

benzo[f]cinnoline has N at positions 4 and 5. The gamma position (3 bonds
from outer N5) is C2: N5→N4→C3→C2. C2-OH converts to benzo[f]cinnolin-2(5H)-one.
The direct alpha C3 (adjacent to N4) does NOT convert.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[f]cinnoline C2 (gamma to N5)
    ("Oc1cnnc2ccc3ccccc3c12",         "benzo[f]cinnolin-2(5H)-one"),
    ("Sc1cnnc2ccc3ccccc3c12",         "benzo[f]cinnolin-2(5H)-thione"),
    # Regression: direct alpha C3 NOT converted
    ("Oc1cc2c(ccc3ccccc32)nn1",       "benzo[f]cinnolin-3-ol"),
    # Regression: parent rings unaffected
    ("c1ccc2c(c1)ccc1nnccc12",        "benzo[f]cinnoline"),
    ("c1ccc2c(c1)nnc1ccccc12",        "benzo[c]cinnoline"),
    # Regression: Phase 768 unchanged
    ("Oc1cnnc2cc3ccccc3cc12",         "benzo[g]cinnolin-5(8H)-one"),
    ("Oc1cnnc2c1ccc1ccccc12",         "benzo[h]cinnolin-7(10H)-one"),
])
def test_phase770_benzo_f_cinnoline_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
