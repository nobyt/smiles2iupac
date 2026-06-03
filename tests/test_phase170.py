"""Phase 170: 単純無機化合物・小分子の保留名追加

IUPAC 2013 P-14.5 等に基づく保留名:
  OO     → hydrogen peroxide
  [H][H] → dihydrogen
  S      → hydrogen sulfide
  O      → water
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("OO",      "hydrogen peroxide"),
    ("[H][H]",  "dihydrogen"),
    ("S",       "hydrogen sulfide"),
    ("O",       "water"),
    # 回帰: 既存の保留名は変わらない
    ("N",       "ammonia"),
    ("N#N",     "dinitrogen"),
    ("O=O",     "dioxygen"),
    ("O=C=O",   "carbon dioxide"),
    ("S=C=S",   "carbon disulfide"),
])
def test_phase170_simple_inorganic_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
