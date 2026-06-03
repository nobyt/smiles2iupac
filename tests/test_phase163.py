"""Phase 163: 炭素・窒素・酸素の単純無機化合物保留名

IUPAC 2013 P-14.5, P-61.5 等に基づく保留名:
  CO2  → carbon dioxide
  CO   → carbon monoxide
  CS2  → carbon disulfide
  N2   → dinitrogen
  O2   → dioxygen
  NH3  → ammonia
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("O=C=O",       "carbon dioxide"),
    ("[C-]#[O+]",   "carbon monoxide"),
    ("S=C=S",       "carbon disulfide"),
    ("N#N",         "dinitrogen"),
    ("O=O",         "dioxygen"),
    ("[NH3]",       "ammonia"),
    # 回帰: 通常の化合物は変わらない
    ("CC(=O)O",     "acetic acid"),
    ("C",           "methane"),
    ("CS(=O)C",     "dimethyl sulfoxide"),
])
def test_phase163_inorganic_special_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
