"""Phase 657: 5,6-/7,8-dihydro variants of quinoxaline/quinazoline/phthalazine/cinnoline (benzene ring partial saturation)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 5,6-dihydroquinoxaline (N1,N4 arom; C5,C6 sp3; C7=C8 sp2)
    ("C1=Cc2nccnc2CC1", "5,6-dihydroquinoxaline"),
    ("Cc1cnc2c(n1)C=CCC2", "2-methyl-5,6-dihydroquinoxaline"),
    ("Cc1cnc2c(n1)CCC=C2", "3-methyl-5,6-dihydroquinoxaline"),
    ("CC1CC=Cc2nccnc21", "5-methyl-5,6-dihydroquinoxaline"),
    ("CC1C=Cc2nccnc2C1", "6-methyl-5,6-dihydroquinoxaline"),
    ("CC1=Cc2nccnc2CC1", "7-methyl-5,6-dihydroquinoxaline"),
    ("CC1=CCCc2nccnc21", "8-methyl-5,6-dihydroquinoxaline"),
    # 5,6-dihydroquinazoline (N1,N3 arom; C5,C6 sp3; C7=C8 sp2)
    ("C1=Cc2ncncc2CC1", "5,6-dihydroquinazoline"),
    ("Cc1ncc2c(n1)C=CCC2", "2-methyl-5,6-dihydroquinazoline"),
    ("Cc1ncnc2c1CCC=C2", "4-methyl-5,6-dihydroquinazoline"),
    ("CC1CC=Cc2ncncc21", "5-methyl-5,6-dihydroquinazoline"),
    ("CC1C=Cc2ncncc2C1", "6-methyl-5,6-dihydroquinazoline"),
    ("CC1=Cc2ncncc2CC1", "7-methyl-5,6-dihydroquinazoline"),
    ("CC1=CCCc2cncnc21", "8-methyl-5,6-dihydroquinazoline"),
    # 7,8-dihydroquinazoline (N1,N3 arom; C7,C8 sp3; C5=C6 sp2)
    ("C1=Cc2cncnc2CC1", "7,8-dihydroquinazoline"),
    ("Cc1ncc2c(n1)CCC=C2", "2-methyl-7,8-dihydroquinazoline"),
    ("Cc1ncnc2c1C=CCC2", "4-methyl-7,8-dihydroquinazoline"),
    ("CC1=CCCc2ncncc21", "5-methyl-7,8-dihydroquinazoline"),
    ("CC1=Cc2cncnc2CC1", "6-methyl-7,8-dihydroquinazoline"),
    ("CC1C=Cc2cncnc2C1", "7-methyl-7,8-dihydroquinazoline"),
    ("CC1CC=Cc2cncnc21", "8-methyl-7,8-dihydroquinazoline"),
    # 5,6-dihydrophthalazine (N2,N3 arom; C5,C6 sp3; C7=C8 sp2)
    ("C1=Cc2cnncc2CC1", "5,6-dihydrophthalazine"),
    ("Cc1nncc2c1C=CCC2", "1-methyl-5,6-dihydrophthalazine"),
    ("Cc1nncc2c1CCC=C2", "4-methyl-5,6-dihydrophthalazine"),
    ("CC1CC=Cc2cnncc21", "5-methyl-5,6-dihydrophthalazine"),
    ("CC1C=Cc2cnncc2C1", "6-methyl-5,6-dihydrophthalazine"),
    ("CC1=Cc2cnncc2CC1", "7-methyl-5,6-dihydrophthalazine"),
    ("CC1=CCCc2cnncc21", "8-methyl-5,6-dihydrophthalazine"),
    # 5,6-dihydrocinnoline (N1,N2 arom; C5,C6 sp3; C7=C8 sp2)
    ("C1=Cc2nnccc2CC1", "5,6-dihydrocinnoline"),
    ("Cc1cc2c(nn1)C=CCC2", "3-methyl-5,6-dihydrocinnoline"),
    ("Cc1cnnc2c1CCC=C2", "4-methyl-5,6-dihydrocinnoline"),
    ("CC1CC=Cc2nnccc21", "5-methyl-5,6-dihydrocinnoline"),
    ("CC1C=Cc2nnccc2C1", "6-methyl-5,6-dihydrocinnoline"),
    ("CC1=Cc2nnccc2CC1", "7-methyl-5,6-dihydrocinnoline"),
    ("CC1=CCCc2ccnnc21", "8-methyl-5,6-dihydrocinnoline"),
    # 7,8-dihydrocinnoline (N1,N2 arom; C7,C8 sp3; C5=C6 sp2)
    ("C1=Cc2ccnnc2CC1", "7,8-dihydrocinnoline"),
    ("Cc1cc2c(nn1)CCC=C2", "3-methyl-7,8-dihydrocinnoline"),
    ("Cc1cnnc2c1C=CCC2", "4-methyl-7,8-dihydrocinnoline"),
    ("CC1=CCCc2nnccc21", "5-methyl-7,8-dihydrocinnoline"),
    ("CC1=Cc2ccnnc2CC1", "6-methyl-7,8-dihydrocinnoline"),
    ("CC1C=Cc2ccnnc2C1", "7-methyl-7,8-dihydrocinnoline"),
    ("CC1CC=Cc2ccnnc21", "8-methyl-7,8-dihydrocinnoline"),
])
def test_phase657(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
