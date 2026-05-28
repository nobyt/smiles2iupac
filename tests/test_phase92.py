"""Phase 92: ジエン/トリエン複数二重結合 (dien/trien 乗数接頭辞)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ジエン + 官能基 (stem に 'a' 保持: IUPAC 2013)
    ("C=CC=CC(=O)C", "hexa-3,5-dien-2-one"),
    ("CC=CC=CC(=O)O", "hexa-2,4-dienoic acid"),  # sorbic acid
    ("C=CC=CC(=O)O", "penta-2,4-dienoic acid"),
    ("C=CC=CCN", "penta-2,4-dien-1-amine"),
    # 単純ジエン (Phase 34 回帰)
    ("C=CC=C", "buta-1,3-diene"),
    ("C=CCC=C", "penta-1,4-diene"),
    ("CC=CC=CC", "hexa-2,4-diene"),
    # 回帰: 単一 ene
    ("CC=CC(=O)C", "pent-3-en-2-one"),
    ("C=CCC(=O)O", "but-3-enoic acid"),
])
def test_phase92_diene_triene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
