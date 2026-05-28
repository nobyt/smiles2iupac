"""Phase 157: アシルアミノ置換基 (acetamido, benzamido 等) (IUPAC 2013 P-66.6.3)

-NH-C(=O)-R 型の N 置換基を {acyl}amido として命名する。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # アシルアミノ基が環の置換基として現れる場合
    ("OC(=O)c1ccc(NC(C)=O)cc1",    "4-acetamidobenzoic acid"),
    ("OC(=O)c1ccc(NC=O)cc1",        "4-formamidobenzoic acid"),
    # 回帰: アミド基が主官能基の場合 (影響を受けないこと)
    ("CC(=O)Nc1ccccc1",            "N-phenylacetamide"),
    ("CC(=O)Nc1ccncc1",            "N-(pyridin-4-yl)acetamide"),
    # 回帰: アミノ基 (NH2)
    ("Nc1ccccc1",                  "aniline"),
    ("Nc1ccncc1",                  "pyridin-4-amine"),
    # 回帰
    ("CC(=O)O",                    "acetic acid"),
])
def test_phase157_acylamino_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
