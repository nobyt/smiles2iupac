"""Phase 101: ヘテロ芳香族 OH/SH → -ol/-thiol サフィックス形式 (IUPAC P-63.1.3, P-63.5.1)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyridinol (-ol suffix)
    ("Oc1ccncc1", "pyridin-4-ol"),
    ("Oc1cccnc1", "pyridin-3-ol"),
    ("Oc1ccccn1", "pyridin-2-ol"),
    # furanols / thiophenols
    ("Oc1ccco1", "furan-2-ol"),
    ("Oc1cccs1", "thiophen-2-ol"),
    # pyrimidinol
    ("Oc1ncccn1", "pyrimidin-2-ol"),
    # 多置換: OH が主官能基で最低ロカント
    ("Cc1ccnc(O)c1", "4-methylpyridin-2-ol"),
    # -thiol suffix (e を省略しない)
    ("Sc1ccncc1", "pyridine-4-thiol"),
    ("Sc1cccs1", "thiophene-2-thiol"),
    # 回帰: benzene 保留名
    ("Oc1ccccc1", "phenol"),
    ("Sc1ccccc1", "benzenethiol"),
    # 回帰: amino suffix も維持
    ("Nc1ccncc1", "pyridin-4-amine"),
])
def test_phase101_heteroaryl_ol_thiol_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
