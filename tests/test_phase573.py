"""Phase 573: Substituted coronene naming.
Coronene: D6h-symmetric, 24 atoms, all 12 CH positions equivalent (locant 1).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61",   "coronene"),
    ("Cc1cc2ccc3ccc4ccc5ccc6ccc1c1c6c5c4c3c21",  "1-methylcoronene"),
])
def test_phase573_coronene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
