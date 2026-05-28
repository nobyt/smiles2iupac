"""Phase 129: ケテン (ethen-1-one, alken-1-one) 命名 (IUPAC 2013 P-66.5.1)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ketene: H₂C=C=O → ethen-1-one
    ("C=C=O", "ethen-1-one"),
    # methylketene: CH₃CH=C=O → prop-1-en-1-one
    ("CC=C=O", "prop-1-en-1-one"),
    # ethylketene
    ("CCC=C=O", "but-1-en-1-one"),
    # 回帰: regular ketones unchanged
    ("CC(=O)C", "propan-2-one"),
    ("CC(=O)CC", "butan-2-one"),
    # 回帰: aldehyde unchanged
    ("CC=O", "ethanal"),
    # 回帰: vinyl aldehyde unchanged (acrolein)
    ("C=CC=O", "prop-2-enal"),
    # 回帰: alkene unchanged
    ("C=C", "ethene"),
])
def test_phase129_ketene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
