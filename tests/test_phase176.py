"""Phase 176: 保留名追加 — セミカルバジド (IUPAC 2013 P-66.5)

  NNC(=O)N → semicarbazide
  NNC(=S)N → thiosemicarbazide
  (N=C=N は Phase 290 で methanediimine に更新)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 新規保留名
    ("N=C=N",       "methanediimine"),
    ("NNC(=O)N",    "semicarbazide"),
    ("NNC(=S)N",    "thiosemicarbazide"),
    # 回帰: 既存保留名は変わらない
    ("NC(=N)N",     "guanidine"),
    ("NC(=O)N",     "urea"),
    ("NC(=S)N",     "thiourea"),
    ("NN",          "hydrazine"),
    ("NNC",         "methylhydrazine"),
])
def test_phase176_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
