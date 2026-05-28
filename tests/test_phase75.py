"""Phase 75: ヒドラジド R-C(=O)-NH-NH₂"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CC(=O)NN", "ethanehydrazide"),
    ("CCC(=O)NN", "propanehydrazide"),
    ("CCCC(=O)NN", "butanehydrazide"),
    ("C(=O)NN", "methanehydrazide"),
    ("CCCCC(=O)NN", "pentanehydrazide"),
])
def test_phase75_hydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
