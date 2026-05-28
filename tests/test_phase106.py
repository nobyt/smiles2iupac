"""Phase 106: 置換ビフェニル命名 (IUPAC P-31.1.3.4, biphenyl retained name)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # パラ置換ビフェニル
    ("Cc1ccc(-c2ccccc2)cc1", "4-methylbiphenyl"),
    ("Clc1ccc(-c2ccccc2)cc1", "4-chlorobiphenyl"),
    ("Fc1ccc(-c2ccccc2)cc1", "4-fluorobiphenyl"),
    ("Brc1ccc(-c2ccccc2)cc1", "4-bromobiphenyl"),
    # オルト・メタ置換ビフェニル
    ("Cc1ccccc1-c1ccccc1", "2-methylbiphenyl"),
    ("Cc1cccc(-c2ccccc2)c1", "3-methylbiphenyl"),
    # SMILES 順序が逆でも同じ (置換環を主環に選択)
    ("c1ccc(-c2ccc(C)cc2)cc1", "4-methylbiphenyl"),
    ("c1ccc(-c2ccc(Cl)cc2)cc1", "4-chlorobiphenyl"),
    # 回帰: 非置換ビフェニル
    ("c1ccc(-c2ccccc2)cc1", "biphenyl"),
    # 回帰: フェニル on ヘテロアリールは変わらず
    ("Cc1ccc(-c2ccncc2)cc1", "4-(4-methylphenyl)pyridine"),
    ("c1ccc(-c2ccncc2)cc1", "4-phenylpyridine"),
])
def test_phase106_substituted_biphenyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
