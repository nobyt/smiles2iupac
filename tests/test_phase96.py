"""Phase 96: ヘテロ芳香族 + アルケン/アルキン置換基 + N-アルキルアミン"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ビニルピリジン (pyridine with ethenyl substituent)
    ("C=Cc1ccncc1", "4-ethenylpyridine"),
    ("C=Cc1cccnc1", "3-ethenylpyridine"),
    ("C=Cc1ccco1", "2-ethenylfuran"),
    # エチニルピリジン
    ("C#Cc1ccncc1", "4-ethynylpyridine"),
    # N-アルキルヘテロアリールアミン
    ("CNc1ccncc1", "N-methylpyridin-4-amine"),
    ("CCNc1ccncc1", "N-ethylpyridin-4-amine"),
    ("CNc1cccnc1", "N-methylpyridin-3-amine"),
    # 回帰: ベンゼン系はそのまま
    ("C=Cc1ccccc1", "styrene"),
    ("C#Cc1ccccc1", "ethynylbenzene"),
    ("CNc1ccccc1", "N-methylaniline"),
])
def test_phase96_heteroaryl_alkenyl_amine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
