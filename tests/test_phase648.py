"""Phase 648: 1,2-dihydroquinoxaline (fix), 1,2-dihydrophthalazine, 1,2-dihydroquinazoline,
3,4-dihydroquinazoline, 1,2,3,4-tetrahydrophthalazine.
Note: sp2 N positions (Hs=0 in parent) are not tested — methylation there breaks C=N,
producing tetrahydro compounds named by a different parent framework.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2-dihydroquinoxaline (renamed from 3,4-; OPSIN-verified; locants 1-3,5-8; not 4=sp2N)
    ("C1=Nc2ccccc2NC1", "1,2-dihydroquinoxaline"),
    ("CN1CC=Nc2ccccc21", "1-methyl-1,2-dihydroquinoxaline"),
    ("CC1C=Nc2ccccc2N1", "2-methyl-1,2-dihydroquinoxaline"),
    ("CC1=Nc2ccccc2NC1", "3-methyl-1,2-dihydroquinoxaline"),
    ("Cc1cccc2c1N=CCN2", "5-methyl-1,2-dihydroquinoxaline"),
    ("Cc1ccc2c(c1)N=CCN2", "6-methyl-1,2-dihydroquinoxaline"),
    ("Cc1ccc2c(c1)NCC=N2", "7-methyl-1,2-dihydroquinoxaline"),
    ("Cc1cccc2c1NCC=N2", "8-methyl-1,2-dihydroquinoxaline"),
    # 1,2-dihydrophthalazine (OPSIN-verified; locants 1,2,4-8; not 3=sp2N)
    ("C1=NNCc2ccccc21", "1,2-dihydrophthalazine"),
    ("CC1NN=Cc2ccccc21", "1-methyl-1,2-dihydrophthalazine"),
    ("CN1Cc2ccccc2C=N1", "2-methyl-1,2-dihydrophthalazine"),
    ("CC1=NNCc2ccccc21", "4-methyl-1,2-dihydrophthalazine"),
    ("Cc1cccc2c1C=NNC2", "5-methyl-1,2-dihydrophthalazine"),
    ("Cc1ccc2c(c1)C=NNC2", "6-methyl-1,2-dihydrophthalazine"),
    ("Cc1ccc2c(c1)CNN=C2", "7-methyl-1,2-dihydrophthalazine"),
    ("Cc1cccc2c1CNN=C2", "8-methyl-1,2-dihydrophthalazine"),
    # 1,2-dihydroquinazoline (OPSIN-verified; locants 1,2,4-8; not 3=sp2N)
    ("C1=NCNc2ccccc21", "1,2-dihydroquinazoline"),
    ("CN1CN=Cc2ccccc21", "1-methyl-1,2-dihydroquinazoline"),
    ("CC1N=Cc2ccccc2N1", "2-methyl-1,2-dihydroquinazoline"),
    ("CC1=NCNc2ccccc21", "4-methyl-1,2-dihydroquinazoline"),
    ("Cc1cccc2c1C=NCN2", "5-methyl-1,2-dihydroquinazoline"),
    ("Cc1ccc2c(c1)C=NCN2", "6-methyl-1,2-dihydroquinazoline"),
    ("Cc1ccc2c(c1)NCN=C2", "7-methyl-1,2-dihydroquinazoline"),
    ("Cc1cccc2c1NCN=C2", "8-methyl-1,2-dihydroquinazoline"),
    # 3,4-dihydroquinazoline (OPSIN-verified; locants 2-8; not 1=sp2N)
    ("C1=Nc2ccccc2CN1", "3,4-dihydroquinazoline"),
    ("CC1=Nc2ccccc2CN1", "2-methyl-3,4-dihydroquinazoline"),
    ("CN1C=Nc2ccccc2C1", "3-methyl-3,4-dihydroquinazoline"),
    ("CC1NC=Nc2ccccc21", "4-methyl-3,4-dihydroquinazoline"),
    ("Cc1cccc2c1CNC=N2", "5-methyl-3,4-dihydroquinazoline"),
    ("Cc1ccc2c(c1)CNC=N2", "6-methyl-3,4-dihydroquinazoline"),
    ("Cc1ccc2c(c1)N=CNC2", "7-methyl-3,4-dihydroquinazoline"),
    ("Cc1cccc2c1N=CNC2", "8-methyl-3,4-dihydroquinazoline"),
    # 1,2,3,4-tetrahydrophthalazine (OPSIN-verified; C2-symmetric: 1=4,2=3,5=8,6=7)
    ("c1ccc2c(c1)CNNC2", "1,2,3,4-tetrahydrophthalazine"),
    ("CC1NNCc2ccccc21", "1-methyl-1,2,3,4-tetrahydrophthalazine"),
    ("CN1Cc2ccccc2CN1", "2-methyl-1,2,3,4-tetrahydrophthalazine"),
    ("Cc1cccc2c1CNNC2", "5-methyl-1,2,3,4-tetrahydrophthalazine"),
    ("Cc1ccc2c(c1)CNNC2", "6-methyl-1,2,3,4-tetrahydrophthalazine"),
])
def test_phase648(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
