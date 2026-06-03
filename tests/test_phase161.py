"""Phase 161: 末端 COOH を持つアルキル置換基の命名

鎖端にカルボン酸基が付く置換基を正確に命名する:
  carboxymethyl   (-CH2-COOH, 炭素数 2, アルキル部 1)
  2-carboxyethyl  (-CH2CH2-COOH, 炭素数 3, アルキル部 2)
  3-carboxypropyl (-CH2CH2CH2-COOH, 炭素数 4, アルキル部 3)

また aggregate_groups の修正により、芳香環 COOH + 鎖 COOH の
dioic_acid 誤集約を防ぐ。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 芳香環 COOH + 鎖 CH2-COOH → 4-(carboxymethyl)benzoic acid
    ("OC(=O)c1ccc(CC(=O)O)cc1",   "4-(carboxymethyl)benzoic acid"),
    # 芳香環 COOH + 鎖 CH2CH2-COOH → 4-(2-carboxyethyl)benzoic acid
    ("OC(=O)c1ccc(CCC(=O)O)cc1",  "4-(2-carboxyethyl)benzoic acid"),
    # 芳香環 COOH + 鎖 CH2CH2CH2-COOH → 4-(3-carboxypropyl)benzoic acid
    ("OC(=O)c1ccc(CCCC(=O)O)cc1", "4-(3-carboxypropyl)benzoic acid"),
    # 回帰: 通常の benzoic acid は変わらない
    ("OC(=O)c1ccccc1",             "benzoic acid"),
    # 回帰: phthalic acid (2 つとも環 COOH, 保留名)
    ("OC(=O)c1ccccc1C(=O)O",       "phthalic acid"),
    # 回帰: 通常のアルキルベンゼンは変わらない
    ("c1ccc(CC)cc1",                "ethylbenzene"),
])
def test_phase161_carboxy_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
