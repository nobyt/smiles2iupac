"""Phase 519: selenoxide / selenone naming (IUPAC 2013 P-65.3.2)

C[Se](=O)C → dimethyl selenoxide; C[Se](=O)(=O)C → dimethyl selenone.
Previously both were mis-detected as selenide.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # selenoxide: C-Se(=O)-C
    ("C[Se](=O)C",       "dimethyl selenoxide"),
    ("CC[Se](=O)CC",     "diethyl selenoxide"),
    ("CC[Se](=O)C",      "ethyl methyl selenoxide"),
    # selenone: C-Se(=O)₂-C
    ("C[Se](=O)(=O)C",   "dimethyl selenone"),
    ("CC[Se](=O)(=O)CC", "diethyl selenone"),
    # regression: plain selenide unaffected
    ("C[Se]C",           "dimethyl selenide"),
    ("CC[SeH]",          "ethaneselenol"),
    # regression: sulfoxide/sulfone unaffected
    ("CS(=O)C",          "dimethyl sulfoxide"),
    ("CS(=O)(=O)C",      "dimethyl sulfone"),
])
def test_phase519_selenoxide_selenone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
