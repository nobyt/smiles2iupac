"""Phase 78: アルケン/アルキン + ニトリル の名前組み立て修正"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ene + nitrile (子音始まり suffix → 'e' を保持; 2炭素は retained name)
    ("C=CC#N", "acrylonitrile"),
    ("C=CCC#N", "but-3-enenitrile"),
    ("C=CCCC#N", "pent-4-enenitrile"),
    # yne + nitrile
    ("C#CCC#N", "but-3-ynenitrile"),
    # 多重結合なし (既存挙動の確認)
    ("CC#N", "acetonitrile"),
    ("CCC#N", "propanenitrile"),
])
def test_phase78_nitrile_elision(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
