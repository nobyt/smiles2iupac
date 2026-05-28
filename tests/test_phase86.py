"""Phase 86: チオ尿素命名 (thiourea, N-methylthiourea 等)
および二級アミン+アルケン (N-methylprop-2-en-1-amine 等)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 二級アミン + ene 鎖
    ("C=CCNC", "N-methylprop-2-en-1-amine"),
    ("C=CCCNC", "N-methylbut-3-en-1-amine"),
    ("CC=CCNC", "N-methylbut-2-en-1-amine"),
    # 回帰: 飽和二級アミン
    ("CNC", "N-methylmethanamine"),
    ("CCNC", "N-methylethanamine"),
    # 回帰: 一次アミン + ene
    ("C=CCN", "prop-2-en-1-amine"),
])
def test_phase86_secondary_amine_with_ene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
