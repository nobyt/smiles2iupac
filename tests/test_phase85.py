"""Phase 85: エステル/酸ハライド/ヒドラジド/チオアミドにアルケン鎖を付加"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # エステル + ene
    ("C=CC(=O)OC", "methyl prop-2-enoate"),
    ("CC=CC(=O)OC", "methyl but-2-enoate"),
    ("C=CC(=O)OCC", "ethyl prop-2-enoate"),
    # 酸ハライド + ene
    ("C=CC(=O)Cl", "prop-2-enoyl chloride"),
    ("C=CC(=O)Br", "prop-2-enoyl bromide"),
    # ヒドラジド + ene
    ("C=CC(=O)NN", "prop-2-enehydrazide"),
    # チオアミド + ene
    ("C=CC(=S)N", "prop-2-enethioamide"),
    # 回帰: 飽和
    ("CC(=O)OC", "methyl acetate"),
    ("CC(=O)Cl", "acetyl chloride"),
    ("CC(=O)NN", "ethanehydrazide"),
    ("CC(=S)N", "ethanethioamide"),
])
def test_phase85_ene_in_naming_paths(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
