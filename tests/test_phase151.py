"""Phase 151: 環状カーボネート・二ヘテロ原子飽和環ケトン (IUPAC 2013)

1,3-dioxolan-2-one (ethylene carbonate), 1,4-dioxane-2-one type lactones,
and base retention names for 1,3-dioxolane, 1,4-dioxane, etc.
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 環状カーボネート (cyclic carbonates)
    ("O=C1OCCO1",         "1,3-dioxolan-2-one"),     # ethylene carbonate
    ("O=C1OCCCO1",        "1,3-dioxan-2-one"),        # trimethylene carbonate
    # 二ヘテロ原子環保留名 (base names)
    ("C1COCO1",           "1,3-dioxolane"),
    ("C1COCCO1",          "1,4-dioxane"),
    ("C1COCOC1",          "1,3-dioxane"),
    ("C1CSCS1",           "1,3-dithiolane"),
    ("C1CSCCO1",          "1,4-oxathiane"),
    # _try_fused_hetero_retained のコア SMILES 正規化修正の回帰
    ("OC(=O)c1ccc2ccccc2n1",  "quinoline-2-carboxylic acid"),
    ("C1COCO1",               "1,3-dioxolane"),
    # 回帰: 既存の官能基
    ("CC(=O)O",               "acetic acid"),
])
def test_phase151_cyclic_carbonates_and_dioxolanes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
