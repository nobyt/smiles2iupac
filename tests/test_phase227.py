"""Phase 227: trialkoxyborane naming (IUPAC 2013 P-68.1.3).

B(OR)3 → tri{alkoxy}borane for symmetric esters.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # symmetric trialkoxyborane
    ("B(OC)(OC)OC",     "trimethoxyborane"),
    ("B(OCC)(OCC)OCC",  "triethoxyborane"),
    # mixed trialkoxyborane
    ("B(OC)(OCC)OCCC",  "(ethoxy)(methoxy)(propoxy)borane"),
    # regression: boronic acid still works
    ("B(O)(O)CC",       "ethylboronic acid"),
])
def test_phase227_trialkoxyborane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
