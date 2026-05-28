"""Phase 159: 縮合ヘテロ芳香族置換基 (indolyl, benzimidazolyl 等) の命名

縮合ヘテロ芳香族がアルキル/酸鎖の置換基として現れる場合、
単一環 (pyrrolyl 等) ではなく縮合環の保留名を使用する。
  e.g. OC(=O)Cc1c[nH]c2ccccc12 → 2-(1H-indol-3-yl)acetic acid
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # インドール置換基
    ("OC(=O)Cc1c[nH]c2ccccc12",        "2-(1H-indol-3-yl)acetic acid"),
    # キノリン置換基
    ("OC(=O)Cc1ccc2ccccc2n1",           "2-(quinolin-2-yl)acetic acid"),
    # ベンゾフラン置換基
    ("OC(=O)Cc1cc2ccccc2o1",            "2-(benzofuran-2-yl)acetic acid"),
    # 単一環ヘテロ芳香族 (回帰)
    ("OC(=O)Cc1ccco1",                  "2-(furan-2-yl)acetic acid"),
    ("OC(=O)Cc1ccncc1",                 "2-(pyridin-4-yl)acetic acid"),
    # 縮合ヘテロ 回帰: 保留名は変わらない
    ("c1ccc2[nH]ccc2c1",                "1H-indole"),
    ("c1ccc2ncccc2c1",                  "quinoline"),
])
def test_phase159_fused_heteroaryl_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
