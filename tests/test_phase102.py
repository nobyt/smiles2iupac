"""Phase 102: 複素芳香族アリール置換基の括弧表記 + 縮合ヘテロ芳香族サフィックス形式"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ヘテロアリール置換基を含むケトン → 括弧付き (IUPAC P-14.5.7)
    ("CC(=O)c1ccncc1", "1-(pyridin-4-yl)ethan-1-one"),
    ("CC(=O)c1cccnc1", "1-(pyridin-3-yl)ethan-1-one"),
    ("CC(=O)c1cccs1", "1-(thiophen-2-yl)ethan-1-one"),
    # 縮合ヘテロ芳香族のアミン/ヒドロキシル → サフィックス形式
    ("Nc1ccc2[nH]ccc2c1", "1H-indol-5-amine"),
    ("Nc1cc2ccccc2[nH]1", "1H-indol-2-amine"),
    ("Oc1ccc2[nH]ccc2c1", "1H-indol-5-ol"),
    ("Nc1ccc2ncccc2c1", "quinolin-6-amine"),
    ("Oc1ccc2ncccc2c1", "quinolin-6-ol"),
    # 回帰: メチル置換縮合ヘテロ芳香族は変わらず
    ("Cc1ccc2[nH]ccc2c1", "5-methyl-1H-indole"),
    ("Cc1ccc2ncccc2c1", "6-methylquinoline"),
    # 回帰: 単純ヘテロ芳香族
    ("Nc1ccncc1", "pyridin-4-amine"),
    ("Oc1ccncc1", "pyridin-4-ol"),
])
def test_phase102_heteroaryl_parens_and_fused_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
