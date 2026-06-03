"""Phase 343: E/Z in thioether (sulfide) parent chain naming (IUPAC 2013).

Unsaturated thioether parent chains now carry E/Z descriptors and double-bond
locants in the substitutive name.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # E/Z in thioether parent chain (alkenyl longer chain)
    ("C/C=C/CSC",                "[(2E)-but-2-en-1-yl] methyl sulfide"),
    (r"C/C=C\CSC",               "[(2Z)-but-2-en-1-yl] methyl sulfide"),
    ("C/C=C/CSCC",               "[(2E)-but-2-en-1-yl] ethyl sulfide"),
    ("CCSC/C=C/C",               "[(2E)-but-2-en-1-yl] ethyl sulfide"),
    # regressions: saturated thioethers unchanged
    ("CSC",                      "dimethyl sulfide"),
    ("CSCC",                     "ethyl methyl sulfide"),
    ("CSCCC",                    "methyl propyl sulfide"),
    ("CC(C)SC",                  "methyl (propan-2-yl) sulfide"),
])
def test_phase343_ez_thioether_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
