"""Phase 76: ジオイック酸 (dioic acid) の置換基ロカント修正"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2-位置換ブタン二酸 (アスパラギン酸は IUPAC 2013 P-12.1 保留名)
    ("NC(CC(=O)O)C(=O)O", "aspartic acid"),
    ("OC(CC(=O)O)C(=O)O", "malic acid"),
    ("ClC(CC(=O)O)C(=O)O", "2-chlorobutanedioic acid"),
    # 非置換は保留名 (IUPAC 2013 P-65.1.1.4)
    ("OC(=O)CC(=O)O", "malonic acid"),
    ("OC(=O)CCC(=O)O", "succinic acid"),
    ("OC(=O)CCCC(=O)O", "glutaric acid"),
    # 2-位置換ペンタン二酸 (グルタミン酸は IUPAC 2013 P-12.1 保留名)
    ("NC(CCC(=O)O)C(=O)O", "glutamic acid"),
])
def test_phase76_dioic_acid_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
