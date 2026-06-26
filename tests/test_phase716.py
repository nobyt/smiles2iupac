"""Phase 716: methyl derivatives of quinoline, isoquinoline, quinolizine,
cinnoline, phthalazine, quinazoline, and quinoxaline.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # quinoline (CH at 2,3,4,5,6,7,8)
    ("c1ccc2ncccc2c1",       "quinoline"),
    ("Cc1ccc2ccccc2n1",      "2-methylquinoline"),
    ("Cc1cnc2ccccc2c1",      "3-methylquinoline"),
    ("Cc1ccnc2ccccc12",      "4-methylquinoline"),
    ("Cc1cccc2ncccc12",      "5-methylquinoline"),
    ("Cc1ccc2ncccc2c1",      "6-methylquinoline"),
    ("Cc1ccc2cccnc2c1",      "7-methylquinoline"),
    ("Cc1cccc2cccnc12",      "8-methylquinoline"),
    # isoquinoline (CH at 1,3,4,5,6,7,8)
    ("c1ccc2cnccc2c1",       "isoquinoline"),
    ("Cc1nccc2ccccc12",      "1-methylisoquinoline"),
    ("Cc1cc2ccccc2cn1",      "3-methylisoquinoline"),
    ("Cc1cncc2ccccc12",      "4-methylisoquinoline"),
    ("Cc1cccc2cnccc12",      "5-methylisoquinoline"),
    ("Cc1ccc2cnccc2c1",      "6-methylisoquinoline"),
    ("Cc1ccc2ccncc2c1",      "7-methylisoquinoline"),
    ("Cc1cccc2ccncc12",      "8-methylisoquinoline"),
    # quinolizine (CH at 1,2,3,4,6,7,8,9)
    ("C1=CCN2C=CC=CC2=C1",   "quinolizine"),
    ("CC1=C2C=CC=CN2CC=C1",  "1-methylquinolizine"),
    ("CC1=CCN2C=CC=CC2=C1",  "2-methylquinolizine"),
    ("CC1=CC=C2C=CC=CN2C1",  "3-methylquinolizine"),
    ("CC1C=CC=C2C=CC=CN21",  "4-methylquinolizine"),
    ("CC1=CC=CC2=CC=CCN12",  "6-methylquinolizine"),
    ("CC1=CN2CC=CC=C2C=C1",  "7-methylquinolizine"),
    ("CC1=CC2=CC=CCN2C=C1",  "8-methylquinolizine"),
    ("CC1=CC=CN2CC=CC=C12",  "9-methylquinolizine"),
    # cinnoline (CH at 3,4,5,6,7,8)
    ("c1ccc2nnccc2c1",       "cinnoline"),
    ("Cc1cc2ccccc2nn1",      "3-methylcinnoline"),
    ("Cc1cnnc2ccccc12",      "4-methylcinnoline"),
    ("Cc1cccc2nnccc12",      "5-methylcinnoline"),
    ("Cc1ccc2nnccc2c1",      "6-methylcinnoline"),
    ("Cc1ccc2ccnnc2c1",      "7-methylcinnoline"),
    ("Cc1cccc2ccnnc12",      "8-methylcinnoline"),
    # phthalazine (CH at 1,5,6)
    ("c1ccc2cnncc2c1",       "phthalazine"),
    ("Cc1nncc2ccccc12",      "1-methylphthalazine"),
    ("Cc1cccc2cnncc12",      "5-methylphthalazine"),
    ("Cc1ccc2cnncc2c1",      "6-methylphthalazine"),
    # quinazoline (CH at 2,4,5,6,7,8)
    ("c1ccc2ncncc2c1",       "quinazoline"),
    ("Cc1ncc2ccccc2n1",      "2-methylquinazoline"),
    ("Cc1ncnc2ccccc12",      "4-methylquinazoline"),
    ("Cc1cccc2ncncc12",      "5-methylquinazoline"),
    ("Cc1ccc2ncncc2c1",      "6-methylquinazoline"),
    ("Cc1ccc2cncnc2c1",      "7-methylquinazoline"),
    ("Cc1cccc2cncnc12",      "8-methylquinazoline"),
    # quinoxaline (CH at 2,5,6)
    ("c1ccc2nccnc2c1",       "quinoxaline"),
    ("Cc1cnc2ccccc2n1",      "2-methylquinoxaline"),
    ("Cc1cccc2nccnc12",      "5-methylquinoxaline"),
    ("Cc1ccc2nccnc2c1",      "6-methylquinoxaline"),
])
def test_phase716(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
