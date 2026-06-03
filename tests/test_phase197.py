"""Phase 197: アルコキシ置換基を持つ炭素の正しい命名

  C(OCC)OCC  → ethoxymethoxyethane  (diethoxymethane の系統名)
  C(OCC)OC   → ethoxymethoxymethane

鎖内炭素に O-C 置換基がある場合、置換基名に O-アルキル部分を含める。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # O-C 置換基を持つ鎖のある化合物
    ("C(OCC)OCC",  "ethoxymethoxyethane"),
    ("C(OCC)OC",   "methoxymethoxyethane"),
    # 回帰: 単純エーテルは変わらない
    ("COC",        "methoxymethane"),
    ("CCOC",       "methoxyethane"),
    ("CCOCC",      "ethoxyethane"),
    # 回帰: アセタールは変わらない
    ("CC(OC)OC",   "1,1-dimethoxyethane"),
    ("CC(OCC)OCC", "1,1-diethoxyethane"),
    ("C(OC)OC",    "dimethoxymethane"),
    # 回帰: アルコキシ置換アルカン
    ("CC(OC)C",    "2-methoxypropane"),
    ("CCOCC",      "ethoxyethane"),
])
def test_phase197_oxy_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
