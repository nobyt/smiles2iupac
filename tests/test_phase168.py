"""Phase 168: アルケン/アルキン ロカント (IUPAC 2013 P-31.1.2.2)

IUPAC 2013 preferred names always cite locants for multiple bonds.
Phase 384 corrected 3-carbon ene/yne to always include the locant:
  C=CC  → prop-1-ene
  CC#C  → prop-1-yne
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 3炭素: IUPAC 2013 はロカントが必要 (Phase 384)
    ("C=CC",    "prop-1-ene"),
    ("CC=C",    "prop-1-ene"),
    ("C#CC",    "prop-1-yne"),
    ("CC#C",    "prop-1-yne"),
    # 4炭素以上: ロカントが必要
    ("C=CCC",   "but-1-ene"),
    ("CC=CC",   "but-2-ene"),
    ("C=CCCC",  "pent-1-ene"),
    ("C#CCCC",  "pent-1-yne"),
    ("CC#CC",   "but-2-yne"),
    # 2炭素: ロカントなし (既存)
    ("C=C",     "ethene"),
    ("C#C",     "ethyne"),
    # 置換プロペン
    ("CC(=C)C", "2-methylprop-1-ene"),
    # 官能基付きプロペン
    ("C=CCO",   "prop-2-en-1-ol"),
])
def test_phase168_propene_propyne(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
