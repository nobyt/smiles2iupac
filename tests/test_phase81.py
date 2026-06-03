"""Phase 81: アルケン+ジオイック酸の名前組み立て修正 (endioic → enedioic)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("OC(=O)C=CC(=O)O", "but-2-enedioic acid"),
    ("OC(=O)C=CCC(=O)O", "pent-2-enedioic acid"),
    ("OC(=O)CC=CC(=O)O", "pent-2-enedioic acid"),
    # 多重結合なし (回帰確認、保留名を使用)
    ("OC(=O)CC(=O)O", "malonic acid"),
    ("OC(=O)CCC(=O)O", "succinic acid"),
])
def test_phase81_enedioic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
