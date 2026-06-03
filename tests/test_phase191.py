"""Phase 191: スルホン酸/スルフィン酸 内部位置ロカントと鎖の正確な選択

  CCCS(=O)(=O)O    → propane-1-sulfonic acid   (not propanesulfonic acid)
  CC(S(=O)(=O)O)C  → propane-2-sulfonic acid   (S on internal C)
  CCC(S(=O)(=O)O)C → butane-2-sulfonic acid    (longest chain through S-C)

IUPAC P-65.3.1: スルホン酸のロカントは常に明示。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # スルホン酸: 末端 S (ロカント = 1)
    ("CCCS(=O)(=O)O",        "propane-1-sulfonic acid"),
    ("CCCCS(=O)(=O)O",       "butane-1-sulfonic acid"),
    # スルホン酸: 内部 S (最長鎖を通る)
    ("CC(S(=O)(=O)O)C",      "propane-2-sulfonic acid"),
    ("CCC(S(=O)(=O)O)C",     "butane-2-sulfonic acid"),
    # スルフィン酸: 同様のロカント
    ("CCCS(=O)O",             "propane-1-sulfinic acid"),
    ("CC(S(=O)(=O)O)C",      "propane-2-sulfonic acid"),
    # スルホンアミド: 3-carbon chain
    ("CCCS(=O)(=O)N",         "propane-1-sulfonamide"),
    # スルホニルクロライド: 3-carbon chain
    ("CCCS(=O)(=O)Cl",        "propane-1-sulfonyl chloride"),
    # 回帰: methane/ethane base → ロカント省略
    ("CS(=O)(=O)O",           "methanesulfonic acid"),
    ("CCS(=O)(=O)O",          "ethanesulfonic acid"),
    ("CS(=O)(=O)N",           "methanesulfonamide"),
])
def test_phase191_sulfonic_acid_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
