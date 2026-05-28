"""Phase 140: 二ヘテロ原子飽和環 (IUPAC 2013 P-31.1.3.4)

1,4-dioxane, 1,3-dioxolane, 1,3-dioxane,
1,3-dithiolane, 1,4-dithiane,
1,4-oxathiane, 1,3-oxathiolane
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # dioxanes
    ("C1COCCO1",  "1,4-dioxane"),
    ("C1OCOCC1",  "1,3-dioxane"),
    # dioxolane (5-membered)
    ("C1COCO1",   "1,3-dioxolane"),
    # dithianes/dithiolanes
    ("C1CSCCS1",  "1,4-dithiane"),
    ("C1CSCS1",   "1,3-dithiolane"),
    # oxathiane/oxathiolane
    ("C1CSCCO1",  "1,4-oxathiane"),
    ("C1CSCO1",   "1,3-oxathiolane"),
    # 回帰: single-heteroatom rings unchanged
    ("C1CCCO1",   "oxolane"),
    ("C1CCCCO1",  "oxane"),
    ("C1CCNCC1",  "piperidine"),
    ("C1COCCN1",  "morpholine"),
])
def test_phase140_dioxane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
