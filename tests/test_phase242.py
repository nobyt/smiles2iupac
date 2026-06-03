"""Phase 242: arsenic acid compounds (IUPAC 2013 P-67.3.3).

  R-As(=O)(OH)2   → {alkyl}arsonic acid
  R-As(OH)2       → {alkyl}arsonous acid
  R2As(=O)(OH)    → di{alkyl}arsinic acid
  R2As-OH         → di{alkyl}arsinous acid
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # arsonic acid: R-As(=O)(OH)2
    ("C[As](=O)(O)O",     "methylarsonic acid"),
    ("CC[As](=O)(O)O",    "ethylarsonic acid"),
    ("CCC[As](=O)(O)O",   "propylarsonic acid"),
    # arsonous acid: R-As(OH)2
    ("C[As](O)O",         "methylarsonous acid"),
    ("CC[As](O)O",        "ethylarsonous acid"),
    # arsinic acid: R2As(=O)(OH)
    ("C[As](=O)(C)O",     "dimethylarsinic acid"),
    ("CC[As](=O)(CC)O",   "diethylarsinic acid"),
    # arsinous acid: R2As-OH
    ("C[As](C)O",         "dimethylarsinous acid"),
    ("CC[As](CC)O",       "diethylarsinous acid"),
    # regression: phosphorus acids unchanged
    ("CP(=O)(O)O",        "methylphosphonic acid"),
    ("CP(O)O",            "methylphosphonous acid"),
])
def test_phase242_arsenic_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
