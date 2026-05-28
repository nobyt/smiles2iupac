"""Phase 87: チオ尿素命名 (thiourea, N-methylthiourea 等)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 無置換チオ尿素
    ("NC(=S)N", "thiourea"),
    # N-単置換
    ("CNC(=S)N", "N-methylthiourea"),
    ("CCNC(=S)N", "N-ethylthiourea"),
    # N,N'-二置換
    ("CNC(=S)NC", "N,N'-dimethylthiourea"),
    ("CCNC(=S)NCC", "N,N'-diethylthiourea"),
    # N,N-同一 N に二置換
    ("CN(C)C(=S)N", "N,N-dimethylthiourea"),
    # 回帰: チオアミド (N-Nなし、C鎖あり → thiourea ではない)
    ("CC(=S)N", "ethanethioamide"),
    # 回帰: チオセミカルバゾン (N-N ありで除外)
    ("CC(=NNC(N)=S)C", "propan-2-one thiosemicarbazone"),
    # 回帰: 尿素 (=O 版)
    ("NC(=O)N", "urea"),
])
def test_phase87_thiourea(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
