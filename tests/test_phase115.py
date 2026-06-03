"""Phase 115: アゾ化合物命名 (IUPAC P-68.3.4)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 対称ジアリール
    ("c1ccc(N=Nc2ccccc2)cc1", "azobenzene"),
    # 対称ジアルキル
    ("CN=NC", "azomethane"),
    ("CCN=NCC", "azoethane"),
    ("CCCN=NCCC", "azopropane"),
    # 回帰: アジド (N=N=N, 3つのN) は対象外
    ("CN=[N+]=[N-]", "azidomethane"),
    # 回帰: ヒドラゾン (C=N-N) は対象外
    ("CC=NN", "ethanal hydrazone"),
    # 回帰: ヒドラジン (N-N 単結合) は対象外
    ("NNC", "methylhydrazine"),
])
def test_phase115_azo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
