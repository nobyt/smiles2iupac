"""Phase 150: ニトレート/ニトライトエステル (nitrate/nitrite esters, IUPAC 2013 P-65.3.1)

Functional class names: {alkyl} nitrate / {alkyl} nitrite.
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Nitrate esters: R-O-[N+](=O)[O-]
    ("CO[N+](=O)[O-]",       "methyl nitrate"),
    ("CCO[N+](=O)[O-]",      "ethyl nitrate"),
    ("CCCO[N+](=O)[O-]",     "propyl nitrate"),
    ("CCCCO[N+](=O)[O-]",    "butyl nitrate"),
    # Nitrite esters: R-O-N=O
    ("CON=O",                "methyl nitrite"),
    ("CCON=O",               "ethyl nitrite"),
    ("CCCON=O",              "propyl nitrite"),
    ("CCCCON=O",             "butyl nitrite"),
    # 回帰
    ("CC(=O)O",              "acetic acid"),
    ("CC[N+](=O)[O-]",       "nitroethane"),
])
def test_phase150_nitrate_nitrite_esters(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
