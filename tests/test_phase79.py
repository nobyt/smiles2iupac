"""Phase 79: N-置換ヒドラゾン命名"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # アルドヒドラゾン N-置換
    ("CC=NNC", "ethanal N-methylhydrazone"),
    ("CCCC=NNC", "butanal N-methylhydrazone"),
    ("CCC=NNC", "propanal N-methylhydrazone"),
    # ケトヒドラゾン N-置換
    ("CC(=NNC)C", "propan-2-one N-methylhydrazone"),
    ("CCC(=NNC)CC", "pentan-3-one N-methylhydrazone"),
    # 未置換 (既存挙動の確認)
    ("CC=NN", "ethanal hydrazone"),
    ("CC(=NN)C", "propan-2-one hydrazone"),
])
def test_phase79_n_substituted_hydrazone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
