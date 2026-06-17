"""Phase 570: Substituted perylene naming.
Perylene: D2h-symmetric, 20 atoms, 3 unique C environments:
A={1,6,7,12}, B={2,5,8,11}, C={3,4,9,10}. IUPAC uses minimum locant.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2cccc3c4cccc5cccc(c(c1)c23)c54",   "perylene"),
    ("Cc1ccc2cccc3c4cccc5cccc(c1c23)c54",    "1-methylperylene"),
    ("Cc1cc2cccc3c4cccc5cccc(c(c1)c23)c54",  "2-methylperylene"),
    ("Cc1ccc2c3cccc4cccc(c5cccc1c52)c43",    "3-methylperylene"),
])
def test_phase570_perylene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
