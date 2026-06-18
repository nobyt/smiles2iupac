"""Phase 618: non-aromatic and partially aromatic rings — methyl-substituted dioxanes,
dithianes, oxathiolanes, indane, tetralin, thiochroman."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2-dioxane (symmetric: 3=6, 4=5)
    ("C1CCOOC1",         "1,2-dioxane"),
    ("CC1CCCOO1",        "3-methyl-1,2-dioxane"),
    ("CC1CCOOC1",        "4-methyl-1,2-dioxane"),
    # 1,2-dioxolane (symmetric: 3=5)
    ("C1COOC1",          "1,2-dioxolane"),
    ("CC1CCOO1",         "3-methyl-1,2-dioxolane"),
    ("CC1COOC1",         "4-methyl-1,2-dioxolane"),
    # 1,2-dithiane (symmetric: 3=6, 4=5)
    ("C1CCSSC1",         "1,2-dithiane"),
    ("CC1CCCSS1",        "3-methyl-1,2-dithiane"),
    ("CC1CCSSC1",        "4-methyl-1,2-dithiane"),
    # 1,2-dithiolane (symmetric: 3=5)
    ("C1CSSC1",          "1,2-dithiolane"),
    ("CC1CCSS1",         "3-methyl-1,2-dithiolane"),
    ("CC1CSSC1",         "4-methyl-1,2-dithiolane"),
    # 1,3-dioxane (C2,C5 distinct; C4=C6 symmetric)
    ("C1COCOC1",         "1,3-dioxane"),
    ("CC1OCCCO1",        "2-methyl-1,3-dioxane"),
    ("CC1CCOCO1",        "4-methyl-1,3-dioxane"),
    ("CC1COCOC1",        "5-methyl-1,3-dioxane"),
    # 1,3-dioxolane (C2 distinct; C4=C5 symmetric)
    ("C1COCO1",          "1,3-dioxolane"),
    ("CC1OCCO1",         "2-methyl-1,3-dioxolane"),
    ("CC1COCO1",         "4-methyl-1,3-dioxolane"),
    # 1,3-dithiolane (C2 distinct; C4=C5 symmetric)
    ("C1CSCS1",          "1,3-dithiolane"),
    ("CC1SCCS1",         "2-methyl-1,3-dithiolane"),
    ("CC1CSCS1",         "4-methyl-1,3-dithiolane"),
    # 1,3-oxathiolane
    ("C1CSCO1",          "1,3-oxathiolane"),
    ("CC1OCCS1",         "2-methyl-1,3-oxathiolane"),
    ("CC1COCS1",         "4-methyl-1,3-oxathiolane"),
    ("CC1CSCO1",         "5-methyl-1,3-oxathiolane"),
    # 1,4-dioxane (all C equivalent → locant 2)
    ("C1COCCO1",         "1,4-dioxane"),
    ("CC1COCCO1",        "2-methyl-1,4-dioxane"),
    # 1,4-dithiane (all C equivalent → locant 2)
    ("C1CSCCS1",         "1,4-dithiane"),
    ("CC1CSCCS1",        "2-methyl-1,4-dithiane"),
    # 1,4-oxathiane
    ("C1CSCCO1",         "1,4-oxathiane"),
    ("CC1COCCS1",        "3-methyl-1,4-oxathiane"),
    ("CC1CSCCO1",        "2-methyl-1,4-oxathiane"),
    # indane (1=3 symmetric, 4=7 symmetric, 5=6 symmetric)
    ("c1ccc2c(c1)CCC2",  "indane"),
    ("CC1CCc2ccccc21",   "1-methylindane"),
    ("CC1Cc2ccccc2C1",   "2-methylindane"),
    ("Cc1cccc2c1CCC2",   "4-methylindane"),
    ("Cc1ccc2c(c1)CCC2", "5-methylindane"),
    # 1,2,3,4-tetrahydronaphthalene (1=4, 2=3, 5=8, 6=7)
    ("c1ccc2c(c1)CCCC2", "1,2,3,4-tetrahydronaphthalene"),
    ("CC1CCCc2ccccc21",  "1-methyl-1,2,3,4-tetrahydronaphthalene"),
    ("CC1CCc2ccccc2C1",  "2-methyl-1,2,3,4-tetrahydronaphthalene"),
    ("Cc1cccc2c1CCCC2",  "5-methyl-1,2,3,4-tetrahydronaphthalene"),
    ("Cc1ccc2c(c1)CCCC2","6-methyl-1,2,3,4-tetrahydronaphthalene"),
    # thiochroman
    ("c1ccc2c(c1)CCCS2", "thiochroman"),
    ("CC1CCc2ccccc2S1",  "2-methylthiochroman"),
    ("CC1CSc2ccccc2C1",  "3-methylthiochroman"),
    ("CC1CCSc2ccccc21",  "4-methylthiochroman"),
    ("Cc1cccc2c1CCCS2",  "5-methylthiochroman"),
    ("Cc1ccc2c(c1)CCCS2","6-methylthiochroman"),
])
def test_phase618_nonaromatic_rings(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
