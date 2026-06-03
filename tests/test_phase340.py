"""Phase 340: E/Z in isocyanate/isothiocyanate substitutive names and carbodiimide (IUPAC 2013).

Unsaturated aliphatic isocyanate/isothiocyanate chains now produce correct
substitutive names with E/Z descriptors. Carbodiimide N-substituents with
E/Z stereo descriptors get proper outer parentheses.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isocyanate E/Z chain
    ("C/C=C/CN=C=O",                 "(2E)-1-isocyanatobut-2-ene"),
    (r"C/C=C\CN=C=O",                "(2Z)-1-isocyanatobut-2-ene"),
    # isothiocyanate E/Z chain
    ("C/C=C/CN=C=S",                 "(2E)-1-isothiocyanatobut-2-ene"),
    # carbodiimide with E/Z N-substituent
    ("C/C=C/CN=C=NCC",               "N-[(2E)-but-2-en-1-yl]-N'-ethylcarbodiimide"),
    ("CCN=C=NC/C=C/C",               "N-[(2E)-but-2-en-1-yl]-N'-ethylcarbodiimide"),
    # regressions: saturated isocyanate/isothiocyanate unchanged
    ("CN=C=O",                       "isocyanatomethane"),
    ("CCN=C=O",                      "isocyanatoethane"),
    ("CC(C)N=C=O",                   "2-isocyanatopropane"),
    ("CCN=C=S",                      "isothiocyanatoethane"),
    ("CN=C=NC",                      "N,N'-dimethylcarbodiimide"),
    ("CCN=C=NCC",                    "N,N'-diethylcarbodiimide"),
])
def test_phase340_ez_isocyanate_carbodiimide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
