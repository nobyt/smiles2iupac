"""Phase 74: ジアシルハライド (diacid halide) 命名"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 対称ジクロライド
    ("ClC(=O)C(=O)Cl", "ethanedioyl dichloride"),
    ("ClC(=O)CC(=O)Cl", "propanedioyl dichloride"),
    ("ClC(=O)CCC(=O)Cl", "butanedioyl dichloride"),
    ("ClC(=O)CCCC(=O)Cl", "pentanedioyl dichloride"),
    # 対称ジブロミド
    ("BrC(=O)C(=O)Br", "ethanedioyl dibromide"),
    ("BrC(=O)CC(=O)Br", "propanedioyl dibromide"),
    # 混合ハライド (アルファベット順)
    ("ClC(=O)C(=O)Br", "ethanedioyl bromide chloride"),
])
def test_phase74_diacid_halide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
