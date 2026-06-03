"""Phase 343: E/Z in thioether (sulfide) parent chain naming (IUPAC 2013).

Unsaturated thioether parent chains now carry E/Z descriptors and double-bond
locants in the substitutive name.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # E/Z in thioether parent chain (alkenyl longer chain)
    ("C/C=C/CSC",                "(2E)-1-(methylsulfanyl)but-2-ene"),
    (r"C/C=C\CSC",               "(2Z)-1-(methylsulfanyl)but-2-ene"),
    ("C/C=C/CSCC",               "(2E)-1-(ethylsulfanyl)but-2-ene"),
    ("CCSC/C=C/C",               "(2E)-1-(ethylsulfanyl)but-2-ene"),
    # regressions: saturated thioethers unchanged
    ("CSC",                      "(methylsulfanyl)methane"),
    ("CSCC",                     "(methylsulfanyl)ethane"),
    ("CSCCC",                    "(methylsulfanyl)propane"),
    ("CC(C)SC",                  "2-(methylsulfanyl)propane"),
])
def test_phase343_ez_thioether_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
