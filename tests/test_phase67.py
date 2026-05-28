"""Phase 67: スルファミド / スルファミン酸命名 (炭素なし特別ケース)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # スルファミド: H2N-S(=O)2-NH2
    ("NS(=O)(=O)N", "sulfamide"),
    # スルファミン酸: H2N-S(=O)2-OH
    ("NS(=O)(=O)O", "sulfamic acid"),
])
def test_phase67_sulfamide_sulfamic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
