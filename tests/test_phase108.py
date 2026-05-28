"""Phase 108: ジアリールケトン (diphenylmethanone) — 1炭素鎖ケトンのロカント省略"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ジフェニルケトン (diphenylmethanone = benzophenone)
    ("O=C(c1ccccc1)c1ccccc1", "diphenylmethanone"),
    # ジアルキルケトンは従来通り (ロカントあり)
    ("CC(=O)C", "propan-2-one"),
    ("CCC(=O)CC", "pentan-3-one"),
    ("CC(=O)CC", "butan-2-one"),
    # フェニルアルキルケトン (親鎖にロカントあり)
    ("O=C(C)c1ccccc1", "acetophenone"),
    ("O=C(CC)c1ccccc1", "1-phenylpropan-1-one"),
    # 回帰: アルデヒドは変わらず
    ("O=Cc1ccccc1", "benzaldehyde"),
    ("CC=O", "ethanal"),
])
def test_phase108_diphenylmethanone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
