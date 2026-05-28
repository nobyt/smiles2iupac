"""Phase 100: ヘテロ芳香族第一級アミン → -amine サフィックス形式 (IUPAC P-62.2.3.1)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ピリジンアミン
    ("Nc1ccncc1", "pyridin-4-amine"),
    ("Nc1cccnc1", "pyridin-3-amine"),
    ("Nc1ccccn1", "pyridin-2-amine"),
    # フラン・チオフェンアミン
    ("Nc1ccco1", "furan-2-amine"),
    ("Nc1cccs1", "thiophen-2-amine"),
    # ピリミジンアミン
    ("Nc1ncccn1", "pyrimidin-2-amine"),
    ("Nc1ccnc(n1)", "pyrimidin-4-amine"),
    # 多置換: アミノが主官能基として最低ロカントを取得
    ("Cc1cncc(N)c1", "5-methylpyridin-3-amine"),
    ("Cc1ccnc(N)c1", "4-methylpyridin-2-amine"),
    # 回帰: ベンゼン → aniline (保留名)
    ("Nc1ccccc1", "aniline"),
    # 回帰: N-置換ヘテロ芳香族アミン (Phase 96) は維持
    ("CNc1ccncc1", "N-methylpyridin-4-amine"),
])
def test_phase100_heteroaryl_primary_amine_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
