"""Phase 628: partially saturated benzo-fused rings — 1,2-dihydronaphthalene,
beta-tetralone, 2,3-dihydro-1,4-benzodioxine, 3,4-dihydro-2H-1,4-benzoxazine,
3,4-dihydro-2H-1,4-benzothiazine, 3,4-dihydroisoquinoline, 3,4-dihydroquinoline,
and 3,4-dihydroquinoxaline, with methyl substituents (IUPAC 2013).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("C1=Cc2ccccc2CC1",   "1,2-dihydronaphthalene"),
    ("O=C1CCc2ccccc2C1",  "3,4-dihydronaphthalen-2(1H)-one"),
    ("c1ccc2c(c1)OCCO2",  "2,3-dihydro-1,4-benzodioxine"),
    ("c1ccc2c(c1)NCCO2",  "3,4-dihydro-2H-1,4-benzoxazine"),
    ("c1ccc2c(c1)NCCS2",  "3,4-dihydro-2H-1,4-benzothiazine"),
    ("C1=NCCc2ccccc21",   "3,4-dihydroisoquinoline"),
    ("C1=Nc2ccccc2CC1",   "3,4-dihydroquinoline"),
    ("C1=Nc2ccccc2NC1",   "3,4-dihydroquinoxaline"),
    # 1,2-dihydronaphthalene: positions 1–8
    ("CC1CC=Cc2ccccc21",  "1-methyl-1,2-dihydronaphthalene"),
    ("CC1C=Cc2ccccc2C1",  "2-methyl-1,2-dihydronaphthalene"),
    ("CC1=Cc2ccccc2CC1",  "3-methyl-1,2-dihydronaphthalene"),
    ("CC1=CCCc2ccccc21",  "4-methyl-1,2-dihydronaphthalene"),
    ("Cc1cccc2c1C=CCC2",  "5-methyl-1,2-dihydronaphthalene"),
    ("Cc1ccc2c(c1)C=CCC2", "6-methyl-1,2-dihydronaphthalene"),
    ("Cc1ccc2c(c1)CCC=C2", "7-methyl-1,2-dihydronaphthalene"),
    ("Cc1cccc2c1CCC=C2",  "8-methyl-1,2-dihydronaphthalene"),
    # 3,4-dihydronaphthalen-2(1H)-one: positions 1, 3–8
    ("CC1C(=O)CCc2ccccc21", "1-methyl-3,4-dihydronaphthalen-2(1H)-one"),
    ("CC1Cc2ccccc2CC1=O",   "3-methyl-3,4-dihydronaphthalen-2(1H)-one"),
    ("CC1CC(=O)Cc2ccccc21", "4-methyl-3,4-dihydronaphthalen-2(1H)-one"),
    ("Cc1cccc2c1CCC(=O)C2", "5-methyl-3,4-dihydronaphthalen-2(1H)-one"),
    ("Cc1ccc2c(c1)CC(=O)CC2", "7-methyl-3,4-dihydronaphthalen-2(1H)-one"),
    # 2,3-dihydro-1,4-benzodioxine: C2-symmetric (2≡3, 5≡8, 6≡7)
    ("CC1COc2ccccc2O1",   "2-methyl-2,3-dihydro-1,4-benzodioxine"),
    ("Cc1cccc2c1OCCO2",   "5-methyl-2,3-dihydro-1,4-benzodioxine"),
    ("Cc1ccc2c(c1)OCCO2", "6-methyl-2,3-dihydro-1,4-benzodioxine"),
    # 3,4-dihydro-2H-1,4-benzoxazine: positions 2, 3, 5–8
    ("CC1CNc2ccccc2O1",   "2-methyl-3,4-dihydro-2H-1,4-benzoxazine"),
    ("CC1COc2ccccc2N1",   "3-methyl-3,4-dihydro-2H-1,4-benzoxazine"),
    ("Cc1cccc2c1NCCO2",   "5-methyl-3,4-dihydro-2H-1,4-benzoxazine"),
    ("Cc1ccc2c(c1)NCCO2", "6-methyl-3,4-dihydro-2H-1,4-benzoxazine"),
    ("Cc1ccc2c(c1)OCCN2", "7-methyl-3,4-dihydro-2H-1,4-benzoxazine"),
    ("Cc1cccc2c1OCCN2",   "8-methyl-3,4-dihydro-2H-1,4-benzoxazine"),
    # 3,4-dihydroisoquinoline: positions 1, 3–8
    ("CC1=NCCc2ccccc21",  "1-methyl-3,4-dihydroisoquinoline"),
    ("CC1Cc2ccccc2C=N1",  "3-methyl-3,4-dihydroisoquinoline"),
    ("CC1CN=Cc2ccccc21",  "4-methyl-3,4-dihydroisoquinoline"),
    ("Cc1cccc2c1CCN=C2",  "5-methyl-3,4-dihydroisoquinoline"),
    ("Cc1ccc2c(c1)C=NCC2", "7-methyl-3,4-dihydroisoquinoline"),
    ("Cc1cccc2c1C=NCC2",  "8-methyl-3,4-dihydroisoquinoline"),
    # 3,4-dihydroquinoline: positions 2–8
    ("CC1=Nc2ccccc2CC1",  "2-methyl-3,4-dihydroquinoline"),
    ("CC1C=Nc2ccccc2C1",  "3-methyl-3,4-dihydroquinoline"),
    ("CC1CC=Nc2ccccc21",  "4-methyl-3,4-dihydroquinoline"),
    ("Cc1cccc2c1CCC=N2",  "5-methyl-3,4-dihydroquinoline"),
    ("Cc1cccc2c1N=CCC2",  "8-methyl-3,4-dihydroquinoline"),
    # 3,4-dihydroquinoxaline: positions 2, 3, 5–8
    ("CC1=Nc2ccccc2NC1",  "2-methyl-3,4-dihydroquinoxaline"),
    ("CC1C=Nc2ccccc2N1",  "3-methyl-3,4-dihydroquinoxaline"),
    ("Cc1cccc2c1NCC=N2",  "5-methyl-3,4-dihydroquinoxaline"),
    ("Cc1ccc2c(c1)NCC=N2", "6-methyl-3,4-dihydroquinoxaline"),
    ("Cc1ccc2c(c1)N=CCN2", "7-methyl-3,4-dihydroquinoxaline"),
    ("Cc1cccc2c1N=CCN2",  "8-methyl-3,4-dihydroquinoxaline"),
])
def test_phase628_partially_saturated_benzo_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
