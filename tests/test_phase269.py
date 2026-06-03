"""Phase 269: substituted monocyclic 5/6-membered diheteroaromatics (IUPAC 2013).

Covers pyrazole, isoxazole, oxazole, thiazole, isothiazole, and the
triazines/tetrazine with methyl substituents at each position.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrazole substituents (lower locant rule; 3/5 tautomeric equivalence)
    ("Cc1ccn[nH]1", "3-methyl-1H-pyrazole"),
    ("Cc1cn[nH]c1", "4-methyl-1H-pyrazole"),
    # isoxazole (O at 1, N at 2)
    ("Cc1ccon1",    "3-methylisoxazole"),
    ("Cc1cnoc1",    "4-methylisoxazole"),
    ("Cc1ccno1",    "5-methylisoxazole"),
    # oxazole (O at 1, N at 3)
    ("Cc1cocn1",    "4-methyloxazole"),
    ("Cc1cnco1",    "5-methyloxazole"),
    # thiazole (S at 1, N at 3)
    ("Cc1cscn1",    "4-methylthiazole"),
    ("Cc1cncs1",    "5-methylthiazole"),
    # isothiazole (S at 1, N at 2)
    ("Cc1cnsc1",    "4-methylisothiazole"),
    ("Cc1ccns1",    "5-methylisothiazole"),
    # 1,2,3-triazine
    ("Cc1ccnnn1",   "4-methyl-1,2,3-triazine"),
    ("Cc1cnnnc1",   "5-methyl-1,2,3-triazine"),
    # 1,2,4-triazine
    ("Cc1cnncn1",   "5-methyl-1,2,4-triazine"),
    # 1,3,5-triazine (symmetric; substituent at C-2)
    ("Cc1ncncn1",   "2-methyl-1,3,5-triazine"),
    # 1,2,4,5-tetrazine (substituent at C-3)
    ("Cc1nncnn1",   "3-methyl-1,2,4,5-tetrazine"),
    # regressions: unsubstituted retained names
    ("c1cn[nH]c1",  "1H-pyrazole"),
    ("c1cnoc1",     "isoxazole"),
    ("c1cocn1",     "oxazole"),
    ("c1cscn1",     "thiazole"),
    ("c1cnsc1",     "isothiazole"),
    ("c1cnnnc1",    "1,2,3-triazine"),
    ("c1ncncn1",    "1,3,5-triazine"),
])
def test_phase269_diheteroaromatic_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
