"""Phase 69: cyanate/thiocyanate functional-class PIN (IUPAC 2013 P-65.3.1).

R-O-C≡N → {R} cyanate; R-S-C≡N → {R} thiocyanate.
Functional-class names are the PINs (same pattern as isocyanate/isothiocyanate).
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # cyanate PIN
    ("COC#N",   "methyl cyanate"),
    ("CCOC#N",  "ethyl cyanate"),
    ("CCCOC#N", "propyl cyanate"),
    # thiocyanate PIN
    ("CSC#N",   "methyl thiocyanate"),
    ("CCSC#N",  "ethyl thiocyanate"),
    ("CCCSC#N", "propyl thiocyanate"),
])
def test_phase69_cyanate_thiocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
