"""Phase 97: ピリミジン置換基のロカント修正 (N-1,3 系)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2-位置換 (N に挟まれた C)
    ("Cc1ncccn1", "2-methylpyrimidine"),
    ("Clc1ncccn1", "2-chloropyrimidine"),
    ("CNc1ncccn1", "N-methylpyrimidin-2-amine"),
    ("OC(=O)c1ncccn1", "pyrimidine-2-carboxylic acid"),
    # 5-位置換 (N に隣接しない C)
    ("Cc1cncnc1", "5-methylpyrimidine"),
    ("Clc1cncnc1", "5-chloropyrimidine"),
    # 4-位置換
    ("Cc1ccnc(n1)", "4-methylpyrimidine"),
    # 回帰: 親環
    ("c1ncccn1", "pyrimidine"),
    ("c1ccncc1", "pyridine"),
    ("Cc1ccncc1", "4-methylpyridine"),
])
def test_phase97_pyrimidine_locants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
