"""Phase 183: 分岐酸エステル・酸ハライドの正しい主鎖選択

  C=C(C)C(=O)OC  → methyl 2-methylprop-2-enoate  (not methyl but-2-enoate)
  CC(C)C(=O)Cl   → 2-methylpropanoyl chloride     (not butanoyl chloride)
  CC(C)C(=O)OC   → methyl 2-methylpropanoate      (not methyl butanoate)

_collect_acid_chain が全炭素収集 (DFS) から最長路探索に修正された。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 分岐メタクリル酸エステル
    ("C=C(C)C(=O)OC",    "methyl 2-methylprop-2-enoate"),
    ("C=C(C)C(=O)OCC",   "ethyl 2-methylprop-2-enoate"),
    # 分岐エステル (飽和)
    ("CC(C)C(=O)OC",     "methyl 2-methylpropanoate"),
    ("CC(C)C(=O)OCC",    "ethyl 2-methylpropanoate"),
    # 分岐酸ハライド
    ("CC(C)C(=O)Cl",     "2-methylpropanoyl chloride"),
    ("CC(C)C(=O)Br",     "2-methylpropanoyl bromide"),
    # 回帰: 直鎖は変わらない
    ("CC(=O)OC",         "methyl acetate"),
    ("CCCC(=O)OC",       "methyl butanoate"),
    ("CC(=O)Cl",         "acetyl chloride"),
])
def test_phase183_branched_acid_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
