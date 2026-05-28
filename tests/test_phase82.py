"""Phase 82: スルフィン酸エステル / スルホン酸エステル命名"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # スルフィン酸エステル C-S(=O)-O-C
    ("CS(=O)OC", "methyl methanesulfinate"),
    ("CS(=O)OCC", "ethyl methanesulfinate"),
    ("CCS(=O)OC", "methyl ethanesulfinate"),
    ("CCCS(=O)OC", "methyl propanesulfinate"),
    # スルホン酸エステル C-S(=O)₂-O-C
    ("CS(=O)(=O)OC", "methyl methanesulfonate"),
    ("CS(=O)(=O)OCC", "ethyl methanesulfonate"),
    ("CCS(=O)(=O)OC", "methyl ethanesulfonate"),
    # 既存スルホキシド/スルホンの回帰確認
    ("CS(=O)C", "dimethyl sulfoxide"),
    ("CS(=O)(=O)C", "dimethyl sulfone"),
])
def test_phase82_sulfinate_sulfonate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
