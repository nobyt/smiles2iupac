"""Phase 112: ジフェニルアミン/トリフェニルアミン保留名 (IUPAC P-62.2.3.2)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ジフェニルアミン保留優先名
    ("c1ccc(Nc2ccccc2)cc1", "diphenylamine"),
    # トリフェニルアミン保留名
    ("c1ccc(N(c2ccccc2)c2ccccc2)cc1", "triphenylamine"),
    # 回帰: N-アルキルアニリン (フェニル置換なし)
    ("c1ccc(NC)cc1", "N-methylaniline"),
    ("c1ccc(NCC)cc1", "N-ethylaniline"),
    ("c1ccc(N(C)C)cc1", "N,N-dimethylaniline"),
    # 回帰: 一級アニリン
    ("Nc1ccccc1", "aniline"),
])
def test_phase112_diphenylamine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
