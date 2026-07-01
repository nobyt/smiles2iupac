"""Phase 76: ジオイック酸 (dioic acid) の置換基ロカント修正"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2-位置換ブタン二酸 → PINs (Phase 735)
    ("NC(CC(=O)O)C(=O)O", "2-aminobutanedioic acid"),
    ("OC(CC(=O)O)C(=O)O", "2-hydroxybutanedioic acid"),
    ("ClC(CC(=O)O)C(=O)O", "2-chlorobutanedioic acid"),
    # 非置換 → systematic PINs (Phase 735)
    ("OC(=O)CC(=O)O", "malonic acid"),
    ("OC(=O)CCC(=O)O", "butanedioic acid"),
    ("OC(=O)CCCC(=O)O", "pentanedioic acid"),
    # 2-位置換ペンタン二酸 → PIN (Phase 735)
    ("NC(CCC(=O)O)C(=O)O", "2-aminopentanedioic acid"),
])
def test_phase76_dioic_acid_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
