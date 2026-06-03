"""Phase 152: イミダゾリジン環 (imidazolidine, imidazolidine-2,4-dione; -dione は e 保持)

5員飽和 N,N 環の保留名と exo C=O ロカント最小化によるジオン命名。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imidazolidine (hydantoin ring without exo C=O)
    ("C1CNCN1",           "imidazolidine"),
    # imidazolidin-2-one (single exo C=O)
    ("O=C1NCCN1",         "imidazolidin-2-one"),
    # imidazolidine-2,4-dione (hydantoin): N at 1,3; C=O at 2,4 (-dione: no elision)
    ("O=C1NC(=O)CN1",     "imidazolidine-2,4-dione"),
    # 回帰: 尿素検出 (環外ウレア) への誤マッチを防ぐ
    ("CNC(=O)N",          "N-methylurea"),
    ("CNC(=O)NC",         "N,N'-dimethylurea"),
    # 回帰: ピペラジン（6員N,N環）
    ("C1CNCCN1",          "piperazine"),
    # 回帰
    ("CC(=O)O",           "acetic acid"),
])
def test_phase152_imidazolidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
