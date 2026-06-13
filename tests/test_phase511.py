"""Phase 511: polyketone tetraone/pentaone suffix support
(IUPAC 2013 P-66.6.2: alkanetetraone, alkanepentaone naming).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CC(=O)CC(=O)CC(=O)CC(=O)C",        "nonane-2,4,6,8-tetraone"),
    ("CC(=O)CC(=O)CC(=O)CC(=O)CC(=O)C",  "undecane-2,4,6,8,10-pentaone"),
    # Lower polyketones still correct
    ("CC(=O)CC(=O)C",    "pentane-2,4-dione"),
    ("CC(=O)CC(=O)CC(=O)C", "heptane-2,4,6-trione"),
])
def test_phase511_higher_polyketone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
