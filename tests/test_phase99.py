"""Phase 99: ヘテロ芳香族多置換体のロカント最小化"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ピリジン多置換体 (N=1 で最小ロカント)
    ("Cc1cc(Cl)ncc1", "2-chloro-4-methylpyridine"),
    ("Clc1cncc(Cl)c1", "3,5-dichloropyridine"),
    ("Cc1cc(C)ncc1", "2,4-dimethylpyridine"),
    ("Clc1cc(Cl)ccn1", "2,4-dichloropyridine"),
    # ピリミジン多置換体
    ("Cc1cc(Cl)nc(n1)", "4-chloro-6-methylpyrimidine"),
    # 単置換は変わらず
    ("Cc1ccncc1", "4-methylpyridine"),
    ("Clc1ccncc1", "4-chloropyridine"),
    ("Cc1cncnc1", "5-methylpyrimidine"),
    ("Cc1ncccn1", "2-methylpyrimidine"),
    # furan 多置換体
    ("Cc1ccc(Cl)o1", "2-chloro-5-methylfuran"),
])
def test_phase99_heteroaryl_multisubstituent_locants(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
