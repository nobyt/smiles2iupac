"""Phase 340: E/Z in isocyanate/isothiocyanate functional-class names and carbodiimide (IUPAC 2013).

Unsaturated aliphatic isocyanate/isothiocyanate chains produce correct
functional-class names with E/Z descriptors on the alkyl group.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isocyanate with E/Z chain
    ("C/C=C/CN=C=O",                 "(2E)-but-2-en-1-yl isocyanate"),
    (r"C/C=C\CN=C=O",                "(2Z)-but-2-en-1-yl isocyanate"),
    # isothiocyanate with E/Z chain
    ("C/C=C/CN=C=S",                 "(2E)-but-2-en-1-yl isothiocyanate"),
    # carbodiimide with E/Z N-substituent
    ("C/C=C/CN=C=NCC",               "N-[(2E)-but-2-en-1-yl]-N'-ethylcarbodiimide"),
    ("CCN=C=NC/C=C/C",               "N-[(2E)-but-2-en-1-yl]-N'-ethylcarbodiimide"),
    # regressions: saturated isocyanate/isothiocyanate
    ("CN=C=O",                       "methyl isocyanate"),
    ("CCN=C=O",                      "ethyl isocyanate"),
    ("CC(C)N=C=O",                   "propan-2-yl isocyanate"),
    ("CCN=C=S",                      "ethyl isothiocyanate"),
    ("CN=C=NC",                      "N,N'-dimethylcarbodiimide"),
    ("CCN=C=NCC",                    "N,N'-diethylcarbodiimide"),
])
def test_phase340_ez_isocyanate_carbodiimide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
