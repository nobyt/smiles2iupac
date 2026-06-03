"""Phase 196: 4員環二ヘテロ原子環・リン含有環 (IUPAC 2013 P-31.1.2)

  C1OC(=O)O1  → 1,3-dioxetan-2-one   (4員環炭酸エステル)
  C1PCC1      → phosphetane           (4員環リン)
  C1NC(=O)O1  → 1,3-oxazetidin-2-one (4員環 O+N)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 4員環 2酸素
    ("C1OC(=O)O1",   "1,3-dioxetan-2-one"),
    ("C1OCO1",       "1,3-dioxetane"),
    # 4員環 O+N
    ("C1NC(=O)O1",   "1,3-oxazetidin-2-one"),
    # 4員環リン (Hantzsch-Widman)
    ("C1PCC1",       "phosphetane"),
    # 3員環リン
    ("C1PC1",        "phosphirane"),
    # 5員環リン
    ("C1CCPC1",      "phospholane"),
    # 6員環リン
    ("C1CCCPC1",     "phosphinane"),
    # 回帰: 既存4員環ヘテロ環は変わらない
    ("C1OCC1",       "oxetane"),
    ("C1NCC1",       "azetidine"),
    ("C1SCC1",       "thietane"),
    ("C1CC(=O)O1",   "oxetan-2-one"),
    ("C1CC(=O)N1",   "azetidin-2-one"),
    # 回帰: 5員環・6員環二ヘテロ原子は変わらない
    ("C1COCO1",      "1,3-dioxolane"),
    ("C1COCCO1",     "1,4-dioxane"),
    ("C1COCOC1",     "1,3-dioxane"),
])
def test_phase196_4membered_dihetero_phosphorus(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
