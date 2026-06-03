"""Phase 246: decahydronaphthalene retained name (IUPAC 2013 P-31.1.3.4).

Decalin (cis or trans, stereochemistry unspecified) → decahydronaphthalene
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # decahydronaphthalene (decalin)
    ("C1CCC2CCCCC2C1",    "decahydronaphthalene"),
    ("C1CCC2CCCCC2CC1",   "decahydronaphthalene"),  # alt SMILES
    # regression: tetralin unchanged
    ("C1CCCc2ccccc12",    "1,2,3,4-tetrahydronaphthalene"),
    # regression: indane unchanged
    ("C1CCc2ccccc21",     "indane"),
    # regression: naphthalene unchanged
    ("c1ccc2ccccc2c1",    "naphthalene"),
    # regression: bicyclo unchanged
    ("C1CC2CCC1C2",       "bicyclo[2.2.1]heptane"),
])
def test_phase246_decahydronaphthalene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
