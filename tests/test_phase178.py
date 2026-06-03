"""Phase 178: 環付きヒドラジド命名 + エン-ジオン命名修正

  O=C(NN)c1ccccc1 → benzohydrazide  (IUPAC 2013 P-65.1.2.4)
  O=C=C=C=O       → propa-1,2-diene-1,3-dione  (e 保持修正)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ベンゾヒドラジド
    ("O=C(NN)c1ccccc1",   "benzohydrazide"),
    # 脂肪族ヒドラジドは変わらない
    ("CC(=O)NN",           "ethanohydrazide"),
    ("CCC(=O)NN",          "propanohydrazide"),
    # カーボンサブオキシド (propa-1,2-diene-1,3-dione)
    ("O=C=C=C=O",          "propa-1,2-diene-1,3-dione"),
    # エン-ジオン (e 保持)
    ("CC(=O)C=CC(=O)C",   "hex-3-ene-2,5-dione"),
    # 回帰: ジオンは変わらない
    ("CC(=O)CC(=O)C",      "pentane-2,4-dione"),
    ("CC(=O)C(=O)C",       "butane-2,3-dione"),
])
def test_phase178_hydrazide_and_diene_dione(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
