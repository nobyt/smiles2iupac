"""Phase 224: sulfinyl chloride (and other sulfinyl halides) naming.

R-S(=O)-X → alkan-N-sulfinyl chloride (IUPAC 2013 P-65.3.1.2).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # sulfinyl chloride
    ("CCS(=O)Cl",   "ethanesulfinyl chloride"),
    ("CS(=O)Cl",    "methanesulfinyl chloride"),
    ("CCCS(=O)Cl",  "propane-1-sulfinyl chloride"),
    # sulfinyl fluoride
    ("CCS(=O)F",    "ethanesulfinyl fluoride"),
    # regression: sulfonyl chloride still works
    ("CCS(=O)(=O)Cl", "ethanesulfonyl chloride"),
])
def test_phase224_sulfinyl_chloride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
