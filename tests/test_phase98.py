"""Phase 98: 置換縮合ヘテロ芳香族 (substituted indole, quinoline, benzimidazole 等)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-indole の置換体
    ("Cc1ccc2[nH]ccc2c1", "5-methyl-1H-indole"),
    ("Cc1cc2ccccc2[nH]1", "2-methyl-1H-indole"),
    ("Cc1c[nH]c2ccccc12", "3-methyl-1H-indole"),
    ("Clc1ccc2[nH]ccc2c1", "5-chloro-1H-indole"),
    # quinoline の置換体
    ("Cc1ccc2ncccc2c1", "6-methylquinoline"),
    ("Clc1ccc2ncccc2c1", "6-chloroquinoline"),
    ("Cc1cccc2ncccc12", "5-methylquinoline"),
    # 1H-benzimidazole の置換体
    ("Cc1ccc2[nH]cnc2c1", "5-methyl-1H-benzimidazole"),
    # benzofuran の置換体
    ("Cc1ccc2occc2c1", "5-methylbenzofuran"),
    # 回帰: 置換なし保留名
    ("c1ccc2[nH]ccc2c1", "1H-indole"),
    ("c1ccc2ncccc2c1", "quinoline"),
    ("c1ccc2[nH]cnc2c1", "1H-benzimidazole"),
    ("c1ccc2occc2c1", "benzofuran"),
    # 回帰: 1H-imidazole の置換
    ("Cc1c[nH]cn1", "4-methyl-1H-imidazole"),
])
def test_phase98_substituted_fused_heteroaromatic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
