"""Phase 166: スルフェン酸 (C-S-OH) の命名

IUPAC 2013 P-65.3 に基づく命名:
  CSO    → methanesulfenic acid
  CCSO   → ethanesulfenic acid
  CCCSO  → propanesulfenic acid
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CSO",          "methanesulfenic acid"),
    ("CCSO",         "ethanesulfenic acid"),
    ("CCCSO",        "propane-1-sulfenic acid"),
    ("CCCCSO",       "butane-1-sulfenic acid"),
    # 回帰: スルフィン酸は変わらない
    ("CS(=O)O",      "methanesulfinic acid"),
    ("CCS(=O)O",     "ethanesulfinic acid"),
    # 回帰: スルホン酸は変わらない
    ("CS(=O)(=O)O",  "methanesulfonic acid"),
    # 回帰: チオールは変わらない
    ("CS",           "methanethiol"),
])
def test_phase166_sulfenic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
