"""Phase 245: organic arsane, organomercury, and inorganic hydrides (IUPAC 2013 P-68).

  R_nAsH_{3-n}   → {alkyl}arsane
  R_nHgH_{2-n}   → di{alkyl}mercury / {halo}({alkyl})mercury
  [AsH3], [GeH4], etc. → retained IUPAC names for inorganic hydrides
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # organic arsanes
    ("C[AsH2]",          "methylarsane"),
    ("CC[AsH2]",         "ethylarsane"),
    ("C[As](C)C",        "trimethylarsane"),
    ("CC[As](CC)CC",     "triethylarsane"),
    # organomercury
    ("C[Hg]C",           "dimethylmercury"),
    ("CC[Hg]CC",         "diethylmercury"),
    ("C[Hg]Cl",          "chloro(methyl)mercury"),
    ("CC[Hg]Br",         "bromo(ethyl)mercury"),
    # inorganic hydrides (no carbon)
    ("[AsH3]",           "arsane"),
    ("[SbH3]",           "stibane"),
    ("[BiH3]",           "bismuthane"),
    ("[GeH4]",           "germane"),
    ("[SnH4]",           "stannane"),
    ("[PbH4]",           "plumbane"),
    # regression: phosphane/stibane unchanged
    ("[PH3]",            "phosphane"),
    ("C[PH2]",           "methylphosphane"),
    ("C[SbH2]",          "methylstibane"),
    # regression: arsenic acids unchanged
    ("C[As](=O)(O)O",    "methylarsonic acid"),
    ("C[As](O)O",        "methylarsonous acid"),
])
def test_phase245_arsane_mercury_hydrides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
