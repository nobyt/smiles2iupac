"""Phase 213: cyclic imine naming (IUPAC 2013 P-66.7.1)

  N=C1CCCCC1  → cyclohexanimine
  N=C1CCCC1   → cyclopentanimine
  N=C1CCC1    → cyclobutanimine

Ring ketones with exocyclic =NH become cycloalkan-1-imines.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Cyclic imines (exocyclic =NH on ring C)
    ("N=C1CCCCC1",  "cyclohexanimine"),
    ("N=C1CCCC1",   "cyclopentanimine"),
    ("N=C1CCC1",    "cyclobutanimine"),
    ("N=C1CC1",     "cyclopropanimine"),
    # regression: acyclic imine unaffected
    ("CC=N",        "ethanimine"),
    ("CC(=N)CC",    "butan-2-imine"),
    # regression: cyclic ketoxime unaffected
    ("ON=C1CCCC1",  "cyclopentanone oxime"),
])
def test_phase213_cyclic_imine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
