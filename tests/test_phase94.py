"""Phase 94: 環状ジアミン + 二級アミン親鎖ロカント修正"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 環状ジアミン
    ("NC1CCCC(N)C1", "cyclohexane-1,3-diamine"),
    ("NC1CCC(N)CC1", "cyclohexane-1,4-diamine"),
    ("NC1CC(N)CC1", "cyclopentane-1,3-diamine"),
    # 二級アミン + 長鎖 (3C以上は locant "-1-" 必要)
    ("CCNCCC", "N-ethylpropan-1-amine"),
    ("CCCNCCC", "N-propylpropan-1-amine"),
    ("CCCCNCC", "N-ethylbutan-1-amine"),
    # 回帰: 短鎖は locant なし
    ("CNC", "N-methylmethanamine"),
    ("CCNC", "N-methylethanamine"),
    ("CCNCC", "N-ethylethanamine"),
    # 回帰: 環状一級アミン
    ("NC1CCCCC1", "cyclohexanamine"),
])
def test_phase94_diamine_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
