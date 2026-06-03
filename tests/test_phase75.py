"""Phase 75: ヒドラジド R-C(=O)-NH-NH₂"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CC(=O)NN", "ethanohydrazide"),
    ("CCC(=O)NN", "propanohydrazide"),
    ("CCCC(=O)NN", "butanohydrazide"),
    ("C(=O)NN", "methanohydrazide"),
    ("CCCCC(=O)NN", "pentanohydrazide"),
])
def test_phase75_hydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
