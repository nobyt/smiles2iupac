"""Phase 108: ジアリールケトン (benzophenone) — 1炭素鎖ケトンのロカント省略"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ジフェニルケトン (benzophenone = retained PIN, IUPAC 2013 P-31.1.3.4)
    ("O=C(c1ccccc1)c1ccccc1", "benzophenone"),
    # アセトン (保留名)
    ("CC(=O)C", "acetone"),
    ("CCC(=O)CC", "pentan-3-one"),
    ("CC(=O)CC", "butan-2-one"),
    # フェニルアルキルケトン (親鎖にロカントあり)
    ("O=C(C)c1ccccc1", "acetophenone"),
    ("O=C(CC)c1ccccc1", "1-phenylpropan-1-one"),
    # 回帰: アルデヒド
    ("O=Cc1ccccc1", "benzaldehyde"),
    ("CC=O", "acetaldehyde"),
])
def test_phase108_diphenylmethanone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
