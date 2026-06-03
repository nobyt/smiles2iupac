"""Phase 93: エンイン/ジイン命名 + ジアルデヒド+ene"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # エンイン (en-yne): IUPAC 2013 P-31.1.6.3 — 二重結合に低いロカントを与える
    ("C=CC#C", "but-1-en-3-yne"),
    ("C#CCC=C", "pent-1-en-4-yne"),
    ("C#CC=C", "but-1-en-3-yne"),
    ("C#CCCC=C", "hex-1-en-5-yne"),
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
