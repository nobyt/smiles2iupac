"""Phase 84: N-置換アミド + アルケン鎖 (N-methylprop-2-enamide 等)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ビニルアミド + N-置換
    ("C=CC(=O)NC", "N-methylprop-2-enamide"),
    ("CC=CC(=O)NC", "N-methylbut-2-enamide"),
    ("C=CC(=O)NCC", "N-ethylprop-2-enamide"),
    ("C=CC(=O)N(C)C", "N,N-dimethylprop-2-enamide"),
    # 回帰: 一次アミドのアルケン (Phase 19)
    ("C=CC(=O)N", "prop-2-enamide"),
    # 回帰: N-置換アミド (二重結合なし)
    ("CC(=O)NC", "N-methylacetamide"),
    ("CC(=O)N(C)C", "N,N-dimethylacetamide"),
])
def test_phase84_n_subst_amide_with_ene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
