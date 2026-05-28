"""Phase 118: Schiff塩基 / N-アルキリデンアリールアミン命名 (IUPAC P-62.3.1)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-エチリデンアニリン
    ("c1ccc(/N=C/C)cc1", "N-ethylideneaniline"),
    ("c1ccc(N=CC)cc1", "N-ethylideneaniline"),
    ("CC=Nc1ccccc1", "N-ethylideneaniline"),
    # N-ベンジリデンアニリン
    ("c1ccc(N=Cc2ccccc2)cc1", "N-benzylideneaniline"),
    # N-メチリデンアニリン (=CH2)
    ("c1ccc(N=C)cc1", "N-methylideneaniline"),
    # 回帰: N-アルキルアニリン (イミンではない)
    ("c1ccc(NC)cc1", "N-methylaniline"),
    ("c1ccc(N)cc1", "aniline"),
    # 回帰: アルキルイミン (N に環なし)
    ("CC(=N)C", "propan-2-imine"),
])
def test_phase118_schiff_base(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
