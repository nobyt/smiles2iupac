"""Phase 149: チオカルボン酸 (thiocarboxylic acids, IUPAC 2013 P-65.1.2.7)

Thioic S-acid (C(=O)SH), thioic O-acid (C(=S)OH), dithioic acid (C(=S)SH).
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Thioic S-acid: R-C(=O)-SH
    ("C(=O)S",           "methanethioic S-acid"),
    ("CC(=O)S",          "ethanethioic S-acid"),
    ("CCC(=O)S",         "propanethioic S-acid"),
    ("CCCC(=O)S",        "butanethioic S-acid"),
    # Thioic O-acid: R-C(=S)-OH
    ("CC(=S)O",          "ethanethioic O-acid"),
    ("CCCC(=S)O",        "butanethioic O-acid"),
    # Dithioic acid: R-C(=S)-SH
    ("CC(=S)S",          "ethanedithioic acid"),
    ("CCCC(=S)S",        "butanedithioic acid"),
    # 回帰: チオアミドは影響なし
    ("CC(=S)N",          "ethanethioamide"),
    ("CC(=O)O",          "acetic acid"),
])
def test_phase149_thioic_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
