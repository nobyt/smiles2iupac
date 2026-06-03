"""Phase 189: ケトン/アルコールの鎖方向決定 - 置換基ロカントによるタイブレーク

  CC(=O)CO → 1-hydroxypropan-2-one  (not 3-hydroxypropan-2-one)
  CC(=O)CBr → 1-bromopropan-2-one   (not 3-bromopropan-2-one)

IUPAC P-14.5: 官能基ロカントが同等の場合、置換基ロカントが最小になる方向を選択。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ケトン中心で置換基ロカントのタイブレーク
    ("CC(=O)CO",     "1-hydroxypropan-2-one"),
    ("CC(=O)CBr",    "1-bromopropan-2-one"),
    ("CC(=O)CCl",    "1-chloropropan-2-one"),
    ("CC(=O)CN",     "1-aminopropan-2-one"),
    # 回帰: 置換基なしケトン
    ("CC(=O)C",      "acetone"),
    ("CC(=O)CC",     "butan-2-one"),
    ("CCC(=O)CC",    "pentan-3-one"),
    # 回帰: 対称ケトンはどちら方向でも同じ
    ("CCC(=O)CCC",   "hexan-3-one"),
])
def test_phase189_ketone_direction_tiebreak(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
