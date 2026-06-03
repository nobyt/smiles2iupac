"""Phase 168: アルケン/アルキン ロカント省略 (IUPAC 2013 P-31.1.2.1)

3炭素鎖の二重結合/三重結合はロカント不要:
  C=CC  → propene  (prop-1-ene ではなく)
  CC#C  → propyne  (prop-1-yne ではなく)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 3炭素: ロカント省略
    ("C=CC",    "propene"),
    ("CC=C",    "propene"),
    ("C#CC",    "propyne"),
    ("CC#C",    "propyne"),
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
