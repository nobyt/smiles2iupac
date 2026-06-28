"""Phase 106: 置換フェニルベンゼン命名 → PIN (IUPAC 2013 P-31.1.3; biphenyl は保留名)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # パラ置換 → X-phenylbenzene
    ("Cc1ccc(-c2ccccc2)cc1", "1-methyl-4-phenylbenzene"),
    ("Clc1ccc(-c2ccccc2)cc1", "1-chloro-4-phenylbenzene"),
    ("Fc1ccc(-c2ccccc2)cc1", "1-fluoro-4-phenylbenzene"),
    ("Brc1ccc(-c2ccccc2)cc1", "1-bromo-4-phenylbenzene"),
    # オルト・メタ置換
    ("Cc1ccccc1-c1ccccc1", "1-methyl-2-phenylbenzene"),
    ("Cc1cccc(-c2ccccc2)c1", "1-methyl-3-phenylbenzene"),
    # SMILES 順序が逆でも同じ
    ("c1ccc(-c2ccc(C)cc2)cc1", "1-methyl-4-phenylbenzene"),
    ("c1ccc(-c2ccc(Cl)cc2)cc1", "1-chloro-4-phenylbenzene"),
    # 非置換ビフェニル → 1,1'-biphenyl
    ("c1ccc(-c2ccccc2)cc1", "1,1'-biphenyl"),
    # 回帰: フェニル on ヘテロアリールは変わらず
    ("Cc1ccc(-c2ccncc2)cc1", "4-(4-methylphenyl)pyridine"),
    ("c1ccc(-c2ccncc2)cc1", "4-phenylpyridine"),
])
def test_phase106_substituted_biphenyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
