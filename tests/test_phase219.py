"""Phase 219: 1,2-dithiolane and 1,2-dithiane (adjacent S-S in ring).

IUPAC 2013: rings with adjacent sulfur atoms at positions 1 and 2.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2-dithiolane: 5-membered ring, S-S adjacent
    ("C1CSSC1",  "1,2-dithiolane"),
    # 1,2-dithiane: 6-membered ring, S-S adjacent
    ("C1CSSCC1", "1,2-dithiane"),
    # regression: 1,3-dithiolane (non-adjacent S) still works
    ("C1CSCS1",  "1,3-dithiolane"),
    # regression: 1,4-dithiane still works
    ("C1CSCCS1", "1,4-dithiane"),
])
def test_phase219_dithiolane_dithiane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
