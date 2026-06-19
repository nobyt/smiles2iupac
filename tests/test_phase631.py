"""Phase 631: tricyclic partially saturated ring systems — 1,2,3,4-tetrahydrophenanthrene,
9,10-dihydrophenanthrene, 1,2-dihydrophenanthrene, 1,2,3,4-tetrahydroanthracene,
9,10-dihydroanthracene, 1,2-dihydroanthracene, and 2,3-dihydro-1H-phenalene,
with methyl substituents at all unique positions (IUPAC 2013).
Also fixes: fluoren-9-one (was 9H-fluoren-9-one).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("c1ccc2c3c(ccc2c1)CCCC3",  "1,2,3,4-tetrahydrophenanthrene"),
    ("c1ccc2c(c1)CCc1ccccc1-2", "9,10-dihydrophenanthrene"),
    ("C1=Cc2c(ccc3ccccc23)CC1", "1,2-dihydrophenanthrene"),
    ("c1ccc2cc3c(cc2c1)CCCC3",  "1,2,3,4-tetrahydroanthracene"),
    ("c1ccc2c(c1)Cc1ccccc1C2",  "9,10-dihydroanthracene"),
    ("C1=Cc2cc3ccccc3cc2CC1",   "1,2-dihydroanthracene"),
    ("c1cc2c3c(cccc3c1)CCC2",   "2,3-dihydro-1H-phenalene"),
    ("O=C1c2ccccc2-c2ccccc21",  "fluoren-9-one"),
    # 1,2,3,4-tetrahydrophenanthrene: positions 1–10
    ("CC1CCCc2c1ccc1ccccc21",   "1-methyl-1,2,3,4-tetrahydrophenanthrene"),
    ("CC1CCc2c(ccc3ccccc23)C1", "2-methyl-1,2,3,4-tetrahydrophenanthrene"),
    ("CC1CCc2ccc3ccccc3c2C1",   "3-methyl-1,2,3,4-tetrahydrophenanthrene"),
    ("CC1CCCc2ccc3ccccc3c21",   "4-methyl-1,2,3,4-tetrahydrophenanthrene"),
    ("Cc1cccc2ccc3c(c12)CCCC3", "5-methyl-1,2,3,4-tetrahydrophenanthrene"),
    ("Cc1ccc2ccc3c(c2c1)CCCC3", "6-methyl-1,2,3,4-tetrahydrophenanthrene"),
    ("Cc1ccc2c3c(ccc2c1)CCCC3", "7-methyl-1,2,3,4-tetrahydrophenanthrene"),
    ("Cc1cccc2c3c(ccc12)CCCC3", "8-methyl-1,2,3,4-tetrahydrophenanthrene"),
    ("Cc1cc2c(c3ccccc13)CCCC2", "9-methyl-1,2,3,4-tetrahydrophenanthrene"),
    ("Cc1cc2ccccc2c2c1CCCC2",   "10-methyl-1,2,3,4-tetrahydrophenanthrene"),
    # 9,10-dihydrophenanthrene: C2-symmetric (1≡6, 2≡5, 3≡4, 9≡10)
    ("Cc1cccc2c1CCc1ccccc1-2",  "1-methyl-9,10-dihydrophenanthrene"),
    ("Cc1ccc2c(c1)CCc1ccccc1-2","2-methyl-9,10-dihydrophenanthrene"),
    ("Cc1ccc2c(c1)-c1ccccc1CC2","3-methyl-9,10-dihydrophenanthrene"),
    ("Cc1cccc2c1-c1ccccc1CC2",  "4-methyl-9,10-dihydrophenanthrene"),
    ("CC1Cc2ccccc2-c2ccccc21",  "9-methyl-9,10-dihydrophenanthrene"),
    # 1,2-dihydrophenanthrene: positions 1–10
    ("CC1CC=Cc2c1ccc1ccccc21",  "1-methyl-1,2-dihydrophenanthrene"),
    ("CC1C=Cc2c(ccc3ccccc23)C1","2-methyl-1,2-dihydrophenanthrene"),
    ("CC1=Cc2c(ccc3ccccc23)CC1","3-methyl-1,2-dihydrophenanthrene"),
    ("CC1=CCCc2ccc3ccccc3c21",  "4-methyl-1,2-dihydrophenanthrene"),
    ("Cc1cccc2ccc3c(c12)C=CCC3","5-methyl-1,2-dihydrophenanthrene"),
    ("Cc1ccc2ccc3c(c2c1)C=CCC3","6-methyl-1,2-dihydrophenanthrene"),
    ("Cc1ccc2c3c(ccc2c1)CCC=C3","7-methyl-1,2-dihydrophenanthrene"),
    ("Cc1cccc2c3c(ccc12)CCC=C3","8-methyl-1,2-dihydrophenanthrene"),
    ("Cc1cc2c(c3ccccc13)C=CCC2","9-methyl-1,2-dihydrophenanthrene"),
    ("Cc1cc2ccccc2c2c1CCC=C2",  "10-methyl-1,2-dihydrophenanthrene"),
    # 1,2,3,4-tetrahydroanthracene: C2-symmetric (1≡4, 2≡3, 5≡9, 6≡... use unique)
    ("CC1CCCc2cc3ccccc3cc21",   "1-methyl-1,2,3,4-tetrahydroanthracene"),
    ("CC1CCc2cc3ccccc3cc2C1",   "2-methyl-1,2,3,4-tetrahydroanthracene"),
    ("Cc1cccc2cc3c(cc12)CCCC3", "5-methyl-1,2,3,4-tetrahydroanthracene"),
    ("Cc1ccc2cc3c(cc2c1)CCCC3", "6-methyl-1,2,3,4-tetrahydroanthracene"),
    ("Cc1c2c(cc3ccccc13)CCCC2", "9-methyl-1,2,3,4-tetrahydroanthracene"),
    # 9,10-dihydroanthracene: C2-symmetric (1≡8, 2≡7, 3≡6, 4≡5, 9≡10)
    ("Cc1cccc2c1Cc1ccccc1C2",   "1-methyl-9,10-dihydroanthracene"),
    ("Cc1ccc2c(c1)Cc1ccccc1C2", "2-methyl-9,10-dihydroanthracene"),
    ("CC1c2ccccc2Cc2ccccc21",   "9-methyl-9,10-dihydroanthracene"),
    # 1,2-dihydroanthracene: positions 1–10
    ("CC1CC=Cc2cc3ccccc3cc21",  "1-methyl-1,2-dihydroanthracene"),
    ("CC1C=Cc2cc3ccccc3cc2C1",  "2-methyl-1,2-dihydroanthracene"),
    ("CC1=Cc2cc3ccccc3cc2CC1",  "3-methyl-1,2-dihydroanthracene"),
    ("CC1=CCCc2cc3ccccc3cc21",  "4-methyl-1,2-dihydroanthracene"),
    ("Cc1cccc2cc3c(cc12)C=CCC3","5-methyl-1,2-dihydroanthracene"),
    ("Cc1ccc2cc3c(cc2c1)C=CCC3","6-methyl-1,2-dihydroanthracene"),
    ("Cc1ccc2cc3c(cc2c1)CCC=C3","7-methyl-1,2-dihydroanthracene"),
    ("Cc1cccc2cc3c(cc12)CCC=C3","8-methyl-1,2-dihydroanthracene"),
    ("Cc1c2c(cc3ccccc13)C=CCC2","9-methyl-1,2-dihydroanthracene"),
    ("Cc1c2c(cc3ccccc13)CCC=C2","10-methyl-1,2-dihydroanthracene"),
    # 2,3-dihydro-1H-phenalene: unique positions 1, 2, 4, 5, 6
    ("CC1CCc2cccc3cccc1c23",    "1-methyl-2,3-dihydro-1H-phenalene"),
    ("CC1Cc2cccc3cccc(c23)C1",  "2-methyl-2,3-dihydro-1H-phenalene"),
    ("Cc1ccc2cccc3c2c1CCC3",    "4-methyl-2,3-dihydro-1H-phenalene"),
    ("Cc1cc2c3c(cccc3c1)CCC2",  "5-methyl-2,3-dihydro-1H-phenalene"),
    ("Cc1ccc2c3c(cccc13)CCC2",  "6-methyl-2,3-dihydro-1H-phenalene"),
])
def test_phase631_tricyclic_partially_saturated(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
