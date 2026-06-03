"""Phase 243: organogermane and organostannane (IUPAC 2013 P-68.1).

  R_nGeH_{4-n}  → {alkyl}germane
  R_nSnH_{4-n}  → {alkyl}stannane
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # organogermanes
    ("C[GeH3]",      "methylgermane"),
    ("CC[GeH3]",     "ethylgermane"),
    ("CC[GeH2]CC",   "diethylgermane"),
    ("C[Ge](C)(C)C", "tetramethylgermane"),
    # organostannanes
    ("C[SnH3]",      "methylstannane"),
    ("CC[SnH3]",     "ethylstannane"),
    ("CC[SnH2]CC",   "diethylstannane"),
    ("C[Sn](C)(C)C", "tetramethylstannane"),
    # regression: silane unchanged
    ("C[SiH3]",      "methylsilane"),
    ("CC[SiH2]CC",   "diethylsilane"),
])
def test_phase243_germane_stannane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
