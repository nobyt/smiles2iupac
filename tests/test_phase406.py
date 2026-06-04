"""Phase 406: Isochromane (isochroman, 3,4-dihydro-1H-2-benzopyran) retained name.

IUPAC 2013 P-31.1.3: isochromane is the retained name for the partially
saturated isobenzopyran ring system where O is at position 2 (not adjacent
to the aromatic ring junction).

Canonical SMILES: c1ccc2c(c1)CCOC2
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # unsubstituted isochromane
    ("C1COCc2ccccc21",                "isochromane"),
    # C-substituted (methyl on benzene ring)
    ("Cc1ccc2c(c1)CCOC2",             "6-methylisochromane"),
    # regression: chromane unchanged
    ("C1CCOc2ccccc21",                "chromane"),
    # regression: indoline unchanged
    ("c1ccc2c(c1)CCN2",               "indoline"),
    # regression: benzene unchanged
    ("c1ccccc1",                       "benzene"),
])
def test_phase406_isochromane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
