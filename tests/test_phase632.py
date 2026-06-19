"""Phase 632: 5,6,7,8-tetrahydro N-heterocycles (aromatic ring intact, benzo ring saturated),
1,2,3,4-tetrahydroquinazoline, quinolizine, and additional polycyclic aromatic hydrocarbons
(benz[a]anthracene, benzo[a/b/c]fluorene, dibenz[a,h]anthracene) — IUPAC 2013.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("c1cnc2c(c1)CCCC2",   "5,6,7,8-tetrahydroquinoline"),
    ("c1cc2c(cn1)CCCC2",   "5,6,7,8-tetrahydroisoquinoline"),
    ("c1cnc2c(n1)CCCC2",   "5,6,7,8-tetrahydroquinoxaline"),
    ("c1ncc2c(n1)CCCC2",   "5,6,7,8-tetrahydroquinazoline"),
    ("c1cc2c(nn1)CCCC2",   "5,6,7,8-tetrahydrocinnoline"),
    ("c1nncc2c1CCCC2",     "5,6,7,8-tetrahydrophthalazine"),
    ("c1ccc2c(c1)CNCN2",   "1,2,3,4-tetrahydroquinazoline"),
    ("C1=CCN2C=CC=CC2=C1", "quinolizine"),
    ("c1ccc2cc3c(ccc4ccccc43)cc2c1",        "benz[a]anthracene"),
    ("c1ccc2c(c1)Cc1c-2ccc2ccccc12",        "benzo[a]fluorene"),
    ("c1ccc2c(c1)Cc1cc3ccccc3cc1-2",        "benzo[b]fluorene"),
    ("c1ccc2c(c1)Cc1ccc3ccccc3c1-2",        "benzo[c]fluorene"),
    ("c1ccc2c(c1)ccc1cc3c(ccc4ccccc43)cc12","dibenz[a,h]anthracene"),
    # 5,6,7,8-tetrahydroquinoline: positions 2–8
    ("Cc1ccc2c(n1)CCCC2",  "2-methyl-5,6,7,8-tetrahydroquinoline"),
    ("Cc1cnc2c(c1)CCCC2",  "3-methyl-5,6,7,8-tetrahydroquinoline"),
    ("Cc1ccnc2c1CCCC2",    "4-methyl-5,6,7,8-tetrahydroquinoline"),
    ("CC1CCCc2ncccc21",    "5-methyl-5,6,7,8-tetrahydroquinoline"),
    ("CC1CCc2ncccc2C1",    "6-methyl-5,6,7,8-tetrahydroquinoline"),
    ("CC1CCc2cccnc2C1",    "7-methyl-5,6,7,8-tetrahydroquinoline"),
    ("CC1CCCc2cccnc21",    "8-methyl-5,6,7,8-tetrahydroquinoline"),
    # 5,6,7,8-tetrahydroisoquinoline: positions 1, 3–8
    ("Cc1nccc2c1CCCC2",    "1-methyl-5,6,7,8-tetrahydroisoquinoline"),
    ("Cc1cc2c(cn1)CCCC2",  "3-methyl-5,6,7,8-tetrahydroisoquinoline"),
    ("Cc1cncc2c1CCCC2",    "4-methyl-5,6,7,8-tetrahydroisoquinoline"),
    ("CC1CCCc2cnccc21",    "5-methyl-5,6,7,8-tetrahydroisoquinoline"),
    ("CC1CCc2cnccc2C1",    "6-methyl-5,6,7,8-tetrahydroisoquinoline"),
    ("CC1CCc2ccncc2C1",    "7-methyl-5,6,7,8-tetrahydroisoquinoline"),
    ("CC1CCCc2ccncc21",    "8-methyl-5,6,7,8-tetrahydroisoquinoline"),
    # 5,6,7,8-tetrahydroquinoxaline: C2-symmetric (2≡3, 5≡8, 6≡7)
    ("Cc1cnc2c(n1)CCCC2",  "2-methyl-5,6,7,8-tetrahydroquinoxaline"),
    ("CC1CCCc2nccnc21",    "5-methyl-5,6,7,8-tetrahydroquinoxaline"),
    ("CC1CCc2nccnc2C1",    "6-methyl-5,6,7,8-tetrahydroquinoxaline"),
    # 5,6,7,8-tetrahydroquinazoline: positions 2, 4–8
    ("Cc1ncc2c(n1)CCCC2",  "2-methyl-5,6,7,8-tetrahydroquinazoline"),
    ("Cc1ncnc2c1CCCC2",    "4-methyl-5,6,7,8-tetrahydroquinazoline"),
    ("CC1CCCc2ncncc21",    "5-methyl-5,6,7,8-tetrahydroquinazoline"),
    ("CC1CCc2ncncc2C1",    "6-methyl-5,6,7,8-tetrahydroquinazoline"),
    ("CC1CCc2cncnc2C1",    "7-methyl-5,6,7,8-tetrahydroquinazoline"),
    ("CC1CCCc2cncnc21",    "8-methyl-5,6,7,8-tetrahydroquinazoline"),
    # 5,6,7,8-tetrahydrocinnoline: positions 3–8
    ("Cc1cc2c(nn1)CCCC2",  "3-methyl-5,6,7,8-tetrahydrocinnoline"),
    ("Cc1cnnc2c1CCCC2",    "4-methyl-5,6,7,8-tetrahydrocinnoline"),
    ("CC1CCCc2nnccc21",    "5-methyl-5,6,7,8-tetrahydrocinnoline"),
    ("CC1CCc2nnccc2C1",    "6-methyl-5,6,7,8-tetrahydrocinnoline"),
    ("CC1CCc2ccnnc2C1",    "7-methyl-5,6,7,8-tetrahydrocinnoline"),
    ("CC1CCCc2ccnnc21",    "8-methyl-5,6,7,8-tetrahydrocinnoline"),
    # 5,6,7,8-tetrahydrophthalazine: C2-symmetric (1≡4, 5≡8, 6≡7)
    ("Cc1nncc2c1CCCC2",    "1-methyl-5,6,7,8-tetrahydrophthalazine"),
    ("CC1CCCc2cnncc21",    "5-methyl-5,6,7,8-tetrahydrophthalazine"),
    ("CC1CCc2cnncc2C1",    "6-methyl-5,6,7,8-tetrahydrophthalazine"),
    # 1,2,3,4-tetrahydroquinazoline: N1-C2-N3-C4-C4a(junc)-C5-C6-C7-C8-C8a(junc)
    ("CN1CNCc2ccccc21",    "1-methyl-1,2,3,4-tetrahydroquinazoline"),
    ("CC1NCc2ccccc2N1",    "2-methyl-1,2,3,4-tetrahydroquinazoline"),
    ("CN1CNc2ccccc2C1",    "3-methyl-1,2,3,4-tetrahydroquinazoline"),
    ("CC1NCNc2ccccc21",    "4-methyl-1,2,3,4-tetrahydroquinazoline"),
    ("Cc1cccc2c1CNCN2",    "5-methyl-1,2,3,4-tetrahydroquinazoline"),
    ("Cc1ccc2c(c1)CNCN2",  "6-methyl-1,2,3,4-tetrahydroquinazoline"),
    ("Cc1ccc2c(c1)NCNC2",  "7-methyl-1,2,3,4-tetrahydroquinazoline"),
    ("Cc1cccc2c1NCNC2",    "8-methyl-1,2,3,4-tetrahydroquinazoline"),
    # quinolizine: positions 1–4, 6–9
    ("CC1=C2C=CC=CN2CC=C1", "1-methylquinolizine"),
    ("CC1=CCN2C=CC=CC2=C1", "2-methylquinolizine"),
    ("CC1=CC=C2C=CC=CN2C1", "3-methylquinolizine"),
    ("CC1C=CC=C2C=CC=CN21", "4-methylquinolizine"),
    ("CC1=CC=CC2=CC=CCN12", "6-methylquinolizine"),
    ("CC1=CN2CC=CC=C2C=C1", "7-methylquinolizine"),
    ("CC1=CC2=CC=CCN2C=C1", "8-methylquinolizine"),
    ("CC1=CC=CN2CC=CC=C12", "9-methylquinolizine"),
])
def test_phase632_tetrahydro_nheterocycles_and_PAH(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
