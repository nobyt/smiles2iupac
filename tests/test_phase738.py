"""Phase 738: ethanohydrazide → acetohydrazide (IUPAC 2013 PIN)

acetic acid (retained PIN) → acetohydrazide, not ethanohydrazide.
N-substituted forms follow the same aceto- stem.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Unsubstituted acetohydrazide
    ("CC(=O)NN",                "acetohydrazide"),
    ("NNC(=O)C",                "acetohydrazide"),
    # N'-substituted
    ("CC(=O)NNC",               "N'-methylacetohydrazide"),
    ("CC(=O)NNCC",              "N'-ethylacetohydrazide"),
    # Regression: C1 and C3+ remain systematic
    ("C(=O)NN",                 "methanohydrazide"),
    ("CCC(=O)NN",               "propanohydrazide"),
    # Regression: benzohydrazide unchanged
    ("O=C(NN)c1ccccc1",         "benzohydrazide"),
])
def test_phase738_acetohydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
