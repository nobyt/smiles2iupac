"""Phase 93: エンイン/ジイン命名 + ジアルデヒド+ene"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # エンイン (en-yne)
    ("C=CC#C", "but-1-en-3-yne"),
    ("C#CCC=C", "pent-4-en-1-yne"),
    ("C#CC=C", "but-3-en-1-yne"),
    ("C#CCCC=C", "hex-5-en-1-yne"),
    # ジイン (diyne)
    ("C#CCCC#C", "hexa-1,5-diyne"),
    ("C#CC#C", "buta-1,3-diyne"),
    # ジアルデヒド + ene ('e' 保持)
    ("O=CC=CC=O", "but-2-enedial"),
    ("O=CC=CC=CC=O", "hexa-2,4-dienedial"),
    # 回帰: 単一アルキン
    ("C#CCC", "but-1-yne"),
    ("CC#CC", "but-2-yne"),
    # 回帰: diene
    ("C=CC=C", "buta-1,3-diene"),
])
def test_phase93_enyne_diyne(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
