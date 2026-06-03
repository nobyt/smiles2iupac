"""Phase 187: ホスフィンオキシド・亜リン酸エステル (IUPAC 2013 P-68)

  CP(=O)(C)C    → trimethylphosphane oxide
  COP(OC)OC     → trimethyl phosphite
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ホスフィンオキシド: R3P=O
    ("CP(=O)(C)C",       "trimethylphosphane oxide"),
    ("O=P(C)(C)C",       "trimethylphosphane oxide"),
    ("CCP(=O)(CC)CC",    "triethylphosphane oxide"),
    ("CP(=O)(C)CC",      "ethyldimethylphosphane oxide"),
    # 亜リン酸トリエステル: (RO)3P
    ("COP(OC)OC",        "trimethyl phosphite"),
    ("CCOP(OCC)OCC",     "triethyl phosphite"),
    # 回帰: ホスファン
    ("CP(C)C",           "trimethylphosphane"),
    ("CCP(CC)CC",        "triethylphosphane"),
    # 回帰: ホスホン酸
    ("CP(=O)(O)O",       "methylphosphonic acid"),
    # 回帰: リン酸エステル
    ("CCOP(=O)(OCC)OCC", "triethyl phosphate"),
])
def test_phase187_phosphine_oxide_phosphite(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
