"""Phase 104: ビアリール化合物の親環選択 — ヘテロ芳香環優先 (IUPAC P-44.1)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # フェニル + ヘテロアリール: ヘテロアリールを親環とする
    ("c1ccc(-c2ccncc2)cc1", "4-phenylpyridine"),
    ("c1ccc(-c2cccnc2)cc1", "3-phenylpyridine"),
    ("c1ccc(-c2ccccn2)cc1", "2-phenylpyridine"),
    ("c1ccc(-c2ccco2)cc1", "2-phenylfuran"),
    ("c1ccc(-c2cccs2)cc1", "2-phenylthiophene"),
    ("c1ccc(-c2cccnc2)cc1", "3-phenylpyridine"),
    # 多置換フェニルピリジン
    ("Cc1ccc(-c2ccncc2)cc1", "4-(4-methylphenyl)pyridine"),
    # 回帰: ビフェニル (2つとも炭素環) は biphenyl
    ("c1ccc(-c2ccccc2)cc1", "biphenyl"),
    # 回帰: 縮合ヘテロ芳香族 (fused) は変わらず
    ("c1ccc2[nH]ccc2c1", "1H-indole"),
    ("c1ccc2ncccc2c1", "quinoline"),
    # 回帰: 置換ヘテロ芳香族
    ("Cc1ccc2[nH]ccc2c1", "5-methyl-1H-indole"),
])
def test_phase104_biaryl_parent_ring_selection(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
