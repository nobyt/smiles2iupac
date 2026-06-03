"""Phase 234: thioketene (C=C=S) → ethene-1-thione (IUPAC 2013 P-65.3.1).

C=C=S has no H on the carbonyl-analog C, so it is detected as thioketone,
not thioaldehyde.  For 2-carbon ene chain, the ene locant is dropped.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C=C=S",     "ethene-1-thione"),
    # longer thioketene analogs
    ("CC=C=S",    "prop-1-ene-1-thione"),
    # regression: normal thioaldehyde unchanged
    ("CCC=S",     "propanethial"),
    ("CC(C)=S",   "propan-2-thione"),
])
def test_phase234_thioketene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
