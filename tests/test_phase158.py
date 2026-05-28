"""Phase 158: N-置換縮合ヘテロ芳香族 (1-methylindole 等) (IUPAC 2013)

N 上に置換基をもつ縮合ヘテロ環: 保留名は 1H- プレフィクスを保持し、
N 位の置換基を "1-{sub}" として命名する。
benzo[b]thiophene への名称修正も含む。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-methyl 縮合ヘテロ芳香族
    ("Cn1ccc2ccccc21",      "1-methyl-1H-indole"),
    ("Cn1cnc2ccccc21",      "1-methyl-1H-benzimidazole"),
    # 回帰: NH 形
    ("c1ccc2[nH]ccc2c1",    "1H-indole"),
    ("c1ccc2[nH]cnc2c1",    "1H-benzimidazole"),
    # benzo[b]thiophene (修正名)
    ("c1ccc2sccc2c1",       "benzo[b]thiophene"),
    # 9H-carbazole (Phase 158 拡張)
    ("c1ccc2[nH]c3ccccc3c2c1",   "9H-carbazole"),
    ("Cn1c2ccccc2c2ccccc21",     "9-methyl-9H-carbazole"),
    # 回帰
    ("c1ccc2ncccc2c1",      "quinoline"),
    ("c1ccc2occc2c1",       "benzofuran"),
    ("CC(=O)O",             "acetic acid"),
])
def test_phase158_n_sub_fused_hetero(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
