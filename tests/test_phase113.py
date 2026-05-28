"""Phase 113: ヒドラジン化合物命名 (IUPAC P-68.3.1)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 親体
    ("NN", "hydrazine"),
    # 一置換アルキル
    ("NNC", "methylhydrazine"),
    ("NNCC", "ethylhydrazine"),
    ("NNCCC", "propylhydrazine"),
    # 一置換アリール
    ("c1ccc(NN)cc1", "phenylhydrazine"),
    # 1,1-二置換
    ("CN(N)C", "1,1-dimethylhydrazine"),
    # 回帰: ヒドラジド (C(=O)-N-N) は対象外 → hydrazide 命名
    ("CC(=O)NN", "ethanehydrazide"),
    # 回帰: アミン
    ("NC", "methanamine"),
    ("Nc1ccccc1", "aniline"),
])
def test_phase113_hydrazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
