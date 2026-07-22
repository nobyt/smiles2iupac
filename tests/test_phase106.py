"""Phase 106/874: 置換ビフェニル命名 → PIN (IUPAC 2013 P-28.2.1)。

二つのベンゼン環からなる環アセンブリの PIN は biphenyl (1,1'-biphenyl) を
親水素化物とする。Phase 106 は当初 "X-phenylbenzene" (置換ベンゼン扱いの
非優先名) を採用していたが、Phase 874 で PIN の "X-...-1,1'-biphenyl" 形式に
修正した (フェニル結合炭素を 1 位に再番号付け)。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # パラ置換 → X-1,1'-biphenyl (Phase 874)
    ("Cc1ccc(-c2ccccc2)cc1", "4-methyl-1,1'-biphenyl"),
    ("Clc1ccc(-c2ccccc2)cc1", "4-chloro-1,1'-biphenyl"),
    ("Fc1ccc(-c2ccccc2)cc1", "4-fluoro-1,1'-biphenyl"),
    ("Brc1ccc(-c2ccccc2)cc1", "4-bromo-1,1'-biphenyl"),
    # オルト・メタ置換
    ("Cc1ccccc1-c1ccccc1", "2-methyl-1,1'-biphenyl"),
    ("Cc1cccc(-c2ccccc2)c1", "3-methyl-1,1'-biphenyl"),
    # SMILES 順序が逆でも同じ
    ("c1ccc(-c2ccc(C)cc2)cc1", "4-methyl-1,1'-biphenyl"),
    ("c1ccc(-c2ccc(Cl)cc2)cc1", "4-chloro-1,1'-biphenyl"),
    # 非置換ビフェニル → 1,1'-biphenyl
    ("c1ccc(-c2ccccc2)cc1", "1,1'-biphenyl"),
    # 回帰: フェニル on ヘテロアリールは変わらず (親=pyridine, biphenyl 経路外)
    ("Cc1ccc(-c2ccncc2)cc1", "4-(4-methylphenyl)pyridine"),
    ("c1ccc(-c2ccncc2)cc1", "4-phenylpyridine"),
])
def test_phase106_substituted_biphenyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
