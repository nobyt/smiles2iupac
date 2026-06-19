"""Phase 627: 2H-indazole, 1,2,3,4-tetrahydroisoquinoline, and 1,2,3,4-tetrahydroquinoxaline
with methyl substituents at all unique positions (IUPAC 2013).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("c1ccc2n[nH]cc2c1",  "2H-indazole"),
    ("c1ccc2c(c1)CCNC2",  "1,2,3,4-tetrahydroisoquinoline"),
    ("c1ccc2c(c1)NCCN2",  "1,2,3,4-tetrahydroquinoxaline"),
    # 2H-indazole: positions 3, 4, 5, 6, 7
    ("Cc1[nH]nc2ccccc12",  "3-methyl-2H-indazole"),
    ("Cc1cccc2n[nH]cc12",  "4-methyl-2H-indazole"),
    ("Cc1ccc2n[nH]cc2c1",  "5-methyl-2H-indazole"),
    ("Cc1ccc2c[nH]nc2c1",  "6-methyl-2H-indazole"),
    ("Cc1cccc2c[nH]nc12",  "7-methyl-2H-indazole"),
    # 1,2,3,4-tetrahydroisoquinoline: positions 1, 3, 4, 5, 6, 7, 8
    ("CC1NCCc2ccccc21",    "1-methyl-1,2,3,4-tetrahydroisoquinoline"),
    ("CC1Cc2ccccc2CN1",    "3-methyl-1,2,3,4-tetrahydroisoquinoline"),
    ("CC1CNCc2ccccc21",    "4-methyl-1,2,3,4-tetrahydroisoquinoline"),
    ("Cc1cccc2c1CCNC2",    "5-methyl-1,2,3,4-tetrahydroisoquinoline"),
    ("Cc1ccc2c(c1)CCNC2",  "6-methyl-1,2,3,4-tetrahydroisoquinoline"),
    ("Cc1ccc2c(c1)CNCC2",  "7-methyl-1,2,3,4-tetrahydroisoquinoline"),
    ("Cc1cccc2c1CNCC2",    "8-methyl-1,2,3,4-tetrahydroisoquinoline"),
    # 1,2,3,4-tetrahydroquinoxaline: C2-symmetric (2≡3, 5≡8, 6≡7)
    ("CC1CNc2ccccc2N1",    "2-methyl-1,2,3,4-tetrahydroquinoxaline"),
    ("Cc1cccc2c1NCCN2",    "5-methyl-1,2,3,4-tetrahydroquinoxaline"),
    ("Cc1ccc2c(c1)NCCN2",  "6-methyl-1,2,3,4-tetrahydroquinoxaline"),
])
def test_phase627_partially_saturated_n_heterocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
