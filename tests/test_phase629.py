"""Phase 629: 1,2-dihydroquinoline and 3,4-dihydroquinoxalin-2(1H)-one
with methyl substituents at all unique positions (IUPAC 2013).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("C1=Cc2ccccc2NC1",   "1,2-dihydroquinoline"),
    ("O=C1CNc2ccccc2N1",  "3,4-dihydroquinoxalin-2(1H)-one"),
    # 1,2-dihydroquinoline: positions 2–8
    ("CC1C=Cc2ccccc2N1",  "2-methyl-1,2-dihydroquinoline"),
    ("CC1=Cc2ccccc2NC1",  "3-methyl-1,2-dihydroquinoline"),
    ("CC1=CCNc2ccccc21",  "4-methyl-1,2-dihydroquinoline"),
    ("Cc1cccc2c1C=CCN2",  "5-methyl-1,2-dihydroquinoline"),
    ("Cc1ccc2c(c1)C=CCN2", "6-methyl-1,2-dihydroquinoline"),
    ("Cc1ccc2c(c1)NCC=C2", "7-methyl-1,2-dihydroquinoline"),
    ("Cc1cccc2c1NCC=C2",  "8-methyl-1,2-dihydroquinoline"),
    # 3,4-dihydroquinoxalin-2(1H)-one: positions 3, 5–8
    ("CC1Nc2ccccc2NC1=O",  "3-methyl-3,4-dihydroquinoxalin-2(1H)-one"),
    ("Cc1cccc2c1NCC(=O)N2", "5-methyl-3,4-dihydroquinoxalin-2(1H)-one"),
    ("Cc1ccc2c(c1)NCC(=O)N2", "6-methyl-3,4-dihydroquinoxalin-2(1H)-one"),
    ("Cc1ccc2c(c1)NC(=O)CN2", "7-methyl-3,4-dihydroquinoxalin-2(1H)-one"),
    ("Cc1cccc2c1NC(=O)CN2", "8-methyl-3,4-dihydroquinoxalin-2(1H)-one"),
])
def test_phase629_dihydroquinoline_quinoxalinone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
