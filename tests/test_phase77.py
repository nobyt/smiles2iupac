"""Phase 77: ペルオキシ酸 R-C(=O)-O-O-H"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C(=O)OO", "methaneperoxoic acid"),
    ("CC(=O)OO", "ethaneperoxoic acid"),
    ("CCC(=O)OO", "propaneperoxoic acid"),
    ("CCCC(=O)OO", "butaneperoxoic acid"),
    ("CCCCC(=O)OO", "pentaneperoxoic acid"),
])
def test_phase77_peroxyacid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
