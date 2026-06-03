"""Phase 177: スルホニルフルオリド命名 (IUPAC 2013 P-65.3.1)

C-S(=O)₂-F → {stem}anesulfonyl fluoride:
  CS(=O)(=O)F           → methanesulfonyl fluoride
  c1ccc(S(=O)(=O)F)cc1  → benzenesulfonyl fluoride
  CS(=O)(=O)Br          → methanesulfonyl bromide
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # スルホニルフルオリド
    ("CS(=O)(=O)F",           "methanesulfonyl fluoride"),
    ("CCS(=O)(=O)F",          "ethanesulfonyl fluoride"),
    ("CCCS(=O)(=O)F",         "propane-1-sulfonyl fluoride"),
    # 芳香族
    ("c1ccc(S(=O)(=O)F)cc1",  "benzenesulfonyl fluoride"),
    # スルホニルブロミド
    ("CS(=O)(=O)Br",          "methanesulfonyl bromide"),
    # 回帰: スルホニルクロライドは変わらない
    ("CS(=O)(=O)Cl",          "methanesulfonyl chloride"),
    ("c1ccc(S(=O)(=O)Cl)cc1", "benzenesulfonyl chloride"),
])
def test_phase177_sulfonyl_halide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
