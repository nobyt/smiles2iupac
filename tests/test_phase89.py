"""Phase 89: スルホン酸/スルホニルクロライド/スルホンアミド+アルケン鎖"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # スルホン酸 + ene
    ("C=CS(=O)(=O)O", "eth-1-enesulfonic acid"),
    ("C=CCS(=O)(=O)O", "prop-2-ene-1-sulfonic acid"),  # allylsulfonic acid
    # スルホニルクロライド + ene
    ("C=CS(=O)(=O)Cl", "eth-1-enesulfonyl chloride"),
    # スルホンアミド + ene
    ("C=CS(=O)(=O)N", "eth-1-enesulfonamide"),
    ("C=CS(=O)(=O)NC", "N-methyleth-1-enesulfonamide"),
    # スルフィン酸 + ene
    ("C=CS(=O)O", "eth-1-enesulfinic acid"),
    # 回帰: 飽和
    ("CS(=O)(=O)O", "methanesulfonic acid"),
    ("CS(=O)(=O)Cl", "methanesulfonyl chloride"),
    ("CS(=O)(=O)N", "methanesulfonamide"),
])
def test_phase89_sulfonate_with_ene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
