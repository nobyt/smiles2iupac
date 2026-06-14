"""Phase 69: cyanate/thiocyanate substitutive PIN (IUPAC 2013 P-65.3.1).

R-O-C≡N → cyanato{alkane}; R-S-C≡N → thiocyanato{alkane}.
Functional class names "methyl cyanate" / "methyl thiocyanate" are retained
acceptable but NOT PINs (contrast: isocyanate R-N=C=O uses functional-class PIN).
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # cyanate PIN
    ("COC#N",   "cyanatomethane"),
    ("CCOC#N",  "cyanatoethane"),
    ("CCCOC#N", "cyanatopropane"),
    # thiocyanate PIN
    ("CSC#N",   "thiocyanatomethane"),
    ("CCSC#N",  "thiocyanatoethane"),
    ("CCCSC#N", "thiocyanatopropane"),
])
def test_phase69_cyanate_thiocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
