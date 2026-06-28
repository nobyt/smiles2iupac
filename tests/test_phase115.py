"""Phase 115: アゾ化合物命名 → diazene PIN (IUPAC 2013 P-68.3.4)
azobenzene/azomethane 等は保留名; PIN は di{R}diazene 形式"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 対称ジアリール → diphenyldiazene
    ("c1ccc(N=Nc2ccccc2)cc1", "diphenyldiazene"),
    # 対称ジアルキル → di{alkyl}diazene
    ("CN=NC", "dimethyldiazene"),
    ("CCN=NCC", "diethyldiazene"),
    ("CCCN=NCCC", "dipropyldiazene"),
    # 回帰: アジド (N=N=N, 3つのN) は対象外
    ("CN=[N+]=[N-]", "azidomethane"),
    # 回帰: ヒドラゾン (C=N-N) は対象外
    ("CC=NN", "ethanal hydrazone"),
    # 回帰: ヒドラジン (N-N 単結合) は対象外
    ("NNC", "methylhydrazine"),
])
def test_phase115_azo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
