"""Phase 638: tricyclic acridine/phenanthridine partially saturated derivatives —
1,2,3,4-tetrahydroacridine, 1,2,3,4-tetrahydrophenanthridine,
9,10-dihydroacridine, and 1,2,3,4,5,6,7,8-octahydroacridine (IUPAC 2013).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("c1ccc2nc3c(cc2c1)CCCC3",  "1,2,3,4-tetrahydroacridine"),
    ("c1ccc2c3c(ncc2c1)CCCC3",  "1,2,3,4-tetrahydrophenanthridine"),
    ("c1ccc2c(c1)Cc1ccccc1N2",  "9,10-dihydroacridine"),
    ("c1c2c(nc3c1CCCC3)CCCC2",  "1,2,3,4,5,6,7,8-octahydroacridine"),
    # 1,2,3,4-tetrahydroacridine: C1-C4(sp3), C4a/C4b/C8a/C9a(junc), N10(None), C5-C8(ar), C9(ar)
    ("CC1CCCc2nc3ccccc3cc21",   "1-methyl-1,2,3,4-tetrahydroacridine"),
    ("CC1CCc2nc3ccccc3cc2C1",   "2-methyl-1,2,3,4-tetrahydroacridine"),
    ("CC1CCc2cc3ccccc3nc2C1",   "3-methyl-1,2,3,4-tetrahydroacridine"),
    ("CC1CCCc2cc3ccccc3nc21",   "4-methyl-1,2,3,4-tetrahydroacridine"),
    ("Cc1cccc2cc3c(nc12)CCCC3", "5-methyl-1,2,3,4-tetrahydroacridine"),
    ("Cc1ccc2cc3c(nc2c1)CCCC3", "6-methyl-1,2,3,4-tetrahydroacridine"),
    ("Cc1ccc2nc3c(cc2c1)CCCC3", "7-methyl-1,2,3,4-tetrahydroacridine"),
    ("Cc1cccc2nc3c(cc12)CCCC3", "8-methyl-1,2,3,4-tetrahydroacridine"),
    ("Cc1c2c(nc3ccccc13)CCCC2", "9-methyl-1,2,3,4-tetrahydroacridine"),
    # 1,2,3,4-tetrahydrophenanthridine: C1-C4(sp3), N5(None), C6-C10(ar)
    ("CC1CCCc2ncc3ccccc3c21",   "1-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("CC1CCc2ncc3ccccc3c2C1",   "2-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("CC1CCc2c(ncc3ccccc23)C1", "3-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("CC1CCCc2c1ncc1ccccc21",   "4-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("Cc1nc2c(c3ccccc13)CCCC2", "6-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("Cc1cccc2c3c(ncc12)CCCC3", "7-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("Cc1ccc2c3c(ncc2c1)CCCC3", "8-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("Cc1ccc2cnc3c(c2c1)CCCC3", "9-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("Cc1cccc2cnc3c(c12)CCCC3", "10-methyl-1,2,3,4-tetrahydrophenanthridine"),
    # 9,10-dihydroacridine: C2-symmetric (1≡5→1, 2≡6→2, 3≡7→3, 4≡8→4); C9(sp3), N10(sp3,NH)
    ("Cc1cccc2c1Cc1ccccc1N2",   "1-methyl-9,10-dihydroacridine"),
    ("Cc1ccc2c(c1)Cc1ccccc1N2", "2-methyl-9,10-dihydroacridine"),
    ("Cc1ccc2c(c1)Nc1ccccc1C2", "3-methyl-9,10-dihydroacridine"),
    ("Cc1cccc2c1Nc1ccccc1C2",   "4-methyl-9,10-dihydroacridine"),
    ("CC1c2ccccc2Nc2ccccc21",   "9-methyl-9,10-dihydroacridine"),
    ("CN1c2ccccc2Cc2ccccc21",   "10-methyl-9,10-dihydroacridine"),
    # 1,2,3,4,5,6,7,8-octahydroacridine: C2-symmetric (1≡5→1, 2≡6→2, 3≡7→3, 4≡8→4); C9(ar,CH)
    ("CC1CCCc2nc3c(cc21)CCCC3", "1-methyl-1,2,3,4,5,6,7,8-octahydroacridine"),
    ("CC1CCc2nc3c(cc2C1)CCCC3", "2-methyl-1,2,3,4,5,6,7,8-octahydroacridine"),
    ("CC1CCc2cc3c(nc2C1)CCCC3", "3-methyl-1,2,3,4,5,6,7,8-octahydroacridine"),
    ("CC1CCCc2cc3c(nc21)CCCC3", "4-methyl-1,2,3,4,5,6,7,8-octahydroacridine"),
    ("Cc1c2c(nc3c1CCCC3)CCCC2", "9-methyl-1,2,3,4,5,6,7,8-octahydroacridine"),
])
def test_phase638_tricyclic_acridines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
