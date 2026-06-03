"""Phase 260: 1,2-dioxane and 1,2-dioxolane — cyclic peroxide ring naming (IUPAC 2013).

  C1CCOOC1 → 1,2-dioxane    (O-O in 6-membered ring, not 'dibutyl peroxide')
  C1COOC1  → 1,2-dioxolane  (O-O in 5-membered ring)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # cyclic O-O compounds (dioxane/dioxolane series)
    ("C1CCOOC1",  "1,2-dioxane"),
    ("C1COOC1",   "1,2-dioxolane"),
    # regression: 1,4-dioxane (non-adjacent O) unchanged
    ("C1COCCO1",  "1,4-dioxane"),
    ("C1COCO1",   "1,3-dioxolane"),
    ("C1COCOC1",  "1,3-dioxane"),
    # regression: open-chain peroxides unchanged
    ("COOC",      "dimethyl peroxide"),
    ("COO",       "methyl hydroperoxide"),
    ("CC(=O)OOC(=O)C", "diethanoyl peroxide"),
])
def test_phase260_dioxane_dioxolane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
