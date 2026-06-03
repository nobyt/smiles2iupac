"""Phase 214: cyclic thioketone naming (IUPAC 2013 P-65.3.1)

  S=C1CCCCC1  → cyclohexanethione
  S=C1CCCC1   → cyclopentanethione

Ring carbons with exocyclic C=S are named as cycloalkanethiones
(loc=1 omitted by convention, analogous to cyclohexanone).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Cyclic thioketones
    ("S=C1CCCCC1",  "cyclohexanethione"),
    ("S=C1CCCC1",   "cyclopentanethione"),
    ("S=C1CCC1",    "cyclobutanethione"),
    ("S=C1CC1",     "cyclopropanethione"),
    # regression: acyclic thioketone unaffected
    ("CC(=S)C",     "propan-2-thione"),
    ("CC(=S)CC",    "butan-2-thione"),
])
def test_phase214_cyclic_thioketone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
