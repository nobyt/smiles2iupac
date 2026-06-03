"""Phase 255: tellurophene, phosphoramides, and boronate esters (IUPAC 2013).

  c1cc[te]c1  → tellurophene (retained heteroaromatic)
  NP(=O)(N)N  → phosphoric triamide
  CB(OC)OC    → dimethyl methylboronate (boronic acid diester)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # tellurophene
    ("c1cc[te]c1",        "tellurophene"),
    # phosphoramides (no carbon, need retained name)
    ("NP(=O)(N)N",        "phosphoric triamide"),
    ("NP(N)N",            "phosphorous triamide"),
    # boronate esters: R-B(OR')2
    ("CB(OC)OC",          "dimethyl methylboronate"),
    ("CCB(OC)OC",         "dimethyl ethylboronate"),
    ("CB(OCC)OCC",        "diethyl methylboronate"),
    # regression: boronic acid unchanged
    ("CB(O)O",            "methylboronic acid"),
    ("CCB(O)O",           "ethylboronic acid"),
    # regression: borate ester unchanged
    ("B(OC)(OC)OC",       "trimethoxyborane"),
    # regression: selenophene unchanged
    ("c1cc[se]c1",        "selenophene"),
])
def test_phase255_tellurophene_boronate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
