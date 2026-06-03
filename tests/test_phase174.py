"""Phase 174: イミン ロカント 1 省略 (IUPAC 2013 P-31.1.2.1)

C=N のみロカント省略 (1C):
  C=N   → methanimine    (not methanimine)
  CC=N  → ethanimine   (pin, locant required)
  CC=NC → N-methylethanimine
3C 以上は区別が必要なためロカントを保持。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1C/2C: ロカント省略
    ("C=N",       "methanimine"),
    ("CC=N",      "ethanimine"),
    ("CC=NC",     "N-methylethanimine"),
    ("CC=NCC",    "N-ethylethanimine"),
    # 3C 以上: ロカント保持
    ("CCC=N",     "propan-1-imine"),
    ("CC(=N)C",   "propan-2-imine"),
    ("CCC=NC",    "N-methylpropan-1-imine"),
    ("CC(=NC)C",  "N-methylpropan-2-imine"),
    # 回帰: 多重結合付きはロカントあり
    ("C=CC=N",    "prop-2-en-1-imine"),
])
def test_phase174_imine_locant_omission(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
