"""Phase 167: リン酸部分エステル (mono/di-ester) の命名修正

IUPAC 2013 P-67.1 に基づく命名:
  CCOP(=O)(O)O   → ethyl dihydrogen phosphate
  CCOP(=O)(OCC)O → diethyl hydrogen phosphate
  COP(=O)(O)O    → methyl dihydrogen phosphate
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # モノエステル
    ("CCOP(=O)(O)O",      "ethyl dihydrogen phosphate"),
    ("COP(=O)(O)O",       "methyl dihydrogen phosphate"),
    # ジエステル
    ("CCOP(=O)(OCC)O",    "diethyl hydrogen phosphate"),
    ("COP(=O)(OC)O",      "dimethyl hydrogen phosphate"),
    # トリエステル (既存)
    ("CCOP(=O)(OCC)OCC",  "triethyl phosphate"),
    ("COP(=O)(OC)OC",     "trimethyl phosphate"),
    # 混合ジエステル
    ("CCOP(=O)(OC)O",     "ethyl methyl hydrogen phosphate"),
    # 回帰: ホスホン酸は変わらない
    ("CP(=O)(O)O",        "methylphosphonic acid"),
])
def test_phase167_phosphate_partial_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
