"""Phase 612: 5-atom monocyclic heteroaromatics — methyl-substituted oxadiazoles,
thiadiazoles, selenophene, tellurophene."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2,5-oxadiazole (symmetric; C3=C4 → locant 3)
    ("c1cnon1",       "1,2,5-oxadiazole"),
    ("Cc1cnon1",      "3-methyl-1,2,5-oxadiazole"),
    # 1,2,5-thiadiazole (symmetric)
    ("c1cnsn1",       "1,2,5-thiadiazole"),
    ("Cc1cnsn1",      "3-methyl-1,2,5-thiadiazole"),
    # 1,2,3-oxadiazole
    ("c1conn1",       "1,2,3-oxadiazole"),
    ("Cc1conn1",      "4-methyl-1,2,3-oxadiazole"),
    ("Cc1cnno1",      "5-methyl-1,2,3-oxadiazole"),
    # 1,2,3-thiadiazole
    ("c1csnn1",       "1,2,3-thiadiazole"),
    ("Cc1csnn1",      "4-methyl-1,2,3-thiadiazole"),
    ("Cc1cnns1",      "5-methyl-1,2,3-thiadiazole"),
    # 1,2,4-oxadiazole
    ("c1ncon1",       "1,2,4-oxadiazole"),
    ("Cc1ncon1",      "3-methyl-1,2,4-oxadiazole"),
    ("Cc1ncno1",      "5-methyl-1,2,4-oxadiazole"),
    # 1,2,4-thiadiazole
    ("c1ncsn1",       "1,2,4-thiadiazole"),
    ("Cc1ncsn1",      "3-methyl-1,2,4-thiadiazole"),
    ("Cc1ncns1",      "5-methyl-1,2,4-thiadiazole"),
    # 1,3,4-oxadiazole (symmetric)
    ("c1nnco1",       "1,3,4-oxadiazole"),
    ("Cc1nnco1",      "2-methyl-1,3,4-oxadiazole"),
    # 1,3,4-thiadiazole (symmetric)
    ("c1nncs1",       "1,3,4-thiadiazole"),
    ("Cc1nncs1",      "2-methyl-1,3,4-thiadiazole"),
    # selenophene
    ("c1cc[se]c1",    "selenophene"),
    ("Cc1ccc[se]1",   "2-methylselenophene"),
    ("Cc1cc[se]c1",   "3-methylselenophene"),
    # tellurophene
    ("c1cc[te]c1",    "tellurophene"),
    ("Cc1ccc[te]1",   "2-methyltellurophene"),
    ("Cc1cc[te]c1",   "3-methyltellurophene"),
])
def test_phase612_monocyclic_heteroaromatics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
