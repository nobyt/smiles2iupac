"""Phase 241: phosphonous acid and dialkylphosphinous acid (IUPAC 2013 P-67.2).

  R-P(OH)2         → {alkyl}phosphonous acid
  R2P-OH           → di{alkyl}phosphinous acid
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phosphonous acid: R-P(OH)2
    ("CP(O)O",       "methylphosphonous acid"),
    ("CCP(O)O",      "ethylphosphonous acid"),
    ("CCCP(O)O",     "propylphosphonous acid"),
    # dialkylphosphinous acid: R2P-OH
    ("CP(C)O",       "dimethylphosphinous acid"),
    ("CCP(CC)O",     "diethylphosphinous acid"),
    # regression: existing phosphinous acid (monoalkyl) unchanged
    ("CP(O)",        "methylphosphinous acid"),
    ("CCP(O)",       "ethylphosphinous acid"),
    # regression: phosphonic acid unchanged
    ("CP(=O)(O)O",   "methylphosphonic acid"),
    # regression: phosphinic acid unchanged
    ("CP(=O)(C)O",   "dimethylphosphinic acid"),
])
def test_phase241_phosphonous_phosphinous(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
