"""Phase 103: ハロゲン化アルキル置換基の命名 + 複合置換基の括弧表記"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ハロゲン化アルキル on ベンゼン
    ("FC(F)(F)c1ccccc1", "(trifluoromethyl)benzene"),
    ("FC(F)c1ccccc1", "(difluoromethyl)benzene"),
    ("ClCc1ccccc1", "(chloromethyl)benzene"),
    ("FC(F)(F)Cc1ccccc1", "(2,2,2-trifluoroethyl)benzene"),
    # ハロゲン化アルキル on ヘテロ芳香環
    ("FC(F)(F)c1ccncc1", "4-(trifluoromethyl)pyridine"),
    ("FC(F)(F)c1cccnc1", "3-(trifluoromethyl)pyridine"),
    ("ClCc1ccncc1", "4-(chloromethyl)pyridine"),
    # ケトン中のヘテロアリール置換基 (Phase 102 との組み合わせ)
    ("CC(=O)c1cccs1", "1-(thiophen-2-yl)ethan-1-one"),
    # 回帰: 単純ハロゲン置換は括弧不要
    ("Clc1ccccc1", "chlorobenzene"),
    ("Clc1ccncc1", "4-chloropyridine"),
    ("Cc1ccncc1", "4-methylpyridine"),
    # 回帰: 主鎖ハロゲン置換は変わらず
    ("CCC(F)(F)F", "1,1,1-trifluoropropane"),
])
def test_phase103_haloalkyl_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
