"""Phase 88: ジエステル + アルケン鎖 (dimethyl but-2-enedioate 等)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ジエステル + ene
    ("COC(=O)C=CC(=O)OC", "dimethyl but-2-enedioate"),
    ("CCOC(=O)C=CC(=O)OCC", "diethyl but-2-enedioate"),
    ("COC(=O)CC=CC(=O)OC", "dimethyl pent-2-enedioate"),
    # 非対称ジエステル + ene
    ("COC(=O)C=CC(=O)OCC", "ethyl methyl but-2-enedioate"),
    # 回帰: 飽和ジエステル (retained names)
    ("COC(=O)CC(=O)OC", "dimethyl malonate"),
    ("CCOC(=O)CC(=O)OCC", "diethyl malonate"),
    ("CCOC(=O)CCC(=O)OCC", "diethyl succinate"),
])
def test_phase88_diester_with_ene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
