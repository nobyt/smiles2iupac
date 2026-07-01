"""Phase 269: substituted monocyclic 5/6-membered diheteroaromatics (IUPAC 2013).

Covers pyrazole, 1,2-oxazole, 1,3-oxazole, 1,3-thiazole, 1,2-thiazole, and the
triazines/tetrazine with methyl substituents at each position.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrazole substituents (lower locant rule; 3/5 tautomeric equivalence)
    ("Cc1cc[nH]n1", "3-methyl-1H-pyrazole"),
    ("Cc1cn[nH]c1", "4-methyl-1H-pyrazole"),
    # 1,2-oxazole (O at 1, N at 2)
    ("Cc1ccon1",    "3-methyl-1,2-oxazole"),
    ("Cc1cnoc1",    "4-methyl-1,2-oxazole"),
    ("Cc1ccno1",    "5-methyl-1,2-oxazole"),
    # 1,3-oxazole (O at 1, N at 3)
    ("Cc1cocn1",    "4-methyl-1,3-oxazole"),
    ("Cc1cnco1",    "5-methyl-1,3-oxazole"),
    # 1,3-thiazole (S at 1, N at 3)
    ("Cc1cscn1",    "4-methyl-1,3-thiazole"),
    ("Cc1cncs1",    "5-methyl-1,3-thiazole"),
    # 1,2-thiazole (S at 1, N at 2)
    ("Cc1cnsc1",    "4-methyl-1,2-thiazole"),
    ("Cc1ccns1",    "5-methyl-1,2-thiazole"),
    # 1,2,3-triazine
    ("Cc1ccnnn1",   "4-methyl-1,2,3-triazine"),
    ("Cc1cnnnc1",   "5-methyl-1,2,3-triazine"),
    # 1,2,4-triazine
    ("Cc1cnncn1",   "5-methyl-1,2,4-triazine"),
    # 1,3,5-triazine (symmetric; substituent at C-2)
    ("Cc1ncncn1",   "2-methyl-1,3,5-triazine"),
    # 1,2,4,5-tetrazine (substituent at C-3)
    ("Cc1nncnn1",   "3-methyl-1,2,4,5-tetrazine"),
    # regressions: unsubstituted PIN names
    ("c1cn[nH]c1",  "1H-pyrazole"),
    ("c1cnoc1",     "1,2-oxazole"),
    ("c1cocn1",     "1,3-oxazole"),
    ("c1cscn1",     "1,3-thiazole"),
    ("c1cnsc1",     "1,2-thiazole"),
    ("c1cnnnc1",    "1,2,3-triazine"),
    ("c1ncncn1",    "1,3,5-triazine"),
])
def test_phase269_diheteroaromatic_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
