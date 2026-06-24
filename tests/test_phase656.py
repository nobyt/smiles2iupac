"""Phase 656: 5,6-/7,8-dihydroquinoline and 5,6-/7,8-dihydroisoquinoline (benzene ring partial saturation)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 7,8-dihydroquinoline (C5=C6 sp2; C7,C8 sp3; N1 arom → no 1-methyl)
    ("C1=Cc2cccnc2CC1", "7,8-dihydroquinoline"),
    ("Cc1ccc2c(n1)CCC=C2", "2-methyl-7,8-dihydroquinoline"),
    ("Cc1cnc2c(c1)C=CCC2", "3-methyl-7,8-dihydroquinoline"),
    ("Cc1ccnc2c1C=CCC2", "4-methyl-7,8-dihydroquinoline"),
    ("CC1=CCCc2ncccc21", "5-methyl-7,8-dihydroquinoline"),
    ("CC1=Cc2cccnc2CC1", "6-methyl-7,8-dihydroquinoline"),
    ("CC1C=Cc2cccnc2C1", "7-methyl-7,8-dihydroquinoline"),
    ("CC1CC=Cc2cccnc21", "8-methyl-7,8-dihydroquinoline"),
    # 5,6-dihydroquinoline (C7=C8 sp2; C5,C6 sp3; N1 arom → no 1-methyl)
    ("C1=Cc2ncccc2CC1", "5,6-dihydroquinoline"),
    ("Cc1ccc2c(n1)C=CCC2", "2-methyl-5,6-dihydroquinoline"),
    ("Cc1cnc2c(c1)CCC=C2", "3-methyl-5,6-dihydroquinoline"),
    ("Cc1ccnc2c1CCC=C2", "4-methyl-5,6-dihydroquinoline"),
    ("CC1CC=Cc2ncccc21", "5-methyl-5,6-dihydroquinoline"),
    ("CC1C=Cc2ncccc2C1", "6-methyl-5,6-dihydroquinoline"),
    ("CC1=Cc2ncccc2CC1", "7-methyl-5,6-dihydroquinoline"),
    ("CC1=CCCc2cccnc21", "8-methyl-5,6-dihydroquinoline"),
    # 5,6-dihydroisoquinoline (C7=C8 sp2; C5,C6 sp3; N2 arom → no 2-methyl)
    ("C1=Cc2cnccc2CC1", "5,6-dihydroisoquinoline"),
    ("Cc1nccc2c1C=CCC2", "1-methyl-5,6-dihydroisoquinoline"),
    ("Cc1cc2c(cn1)C=CCC2", "3-methyl-5,6-dihydroisoquinoline"),
    ("Cc1cncc2c1CCC=C2", "4-methyl-5,6-dihydroisoquinoline"),
    ("CC1CC=Cc2cnccc21", "5-methyl-5,6-dihydroisoquinoline"),
    ("CC1C=Cc2cnccc2C1", "6-methyl-5,6-dihydroisoquinoline"),
    ("CC1=Cc2cnccc2CC1", "7-methyl-5,6-dihydroisoquinoline"),
    ("CC1=CCCc2ccncc21", "8-methyl-5,6-dihydroisoquinoline"),
    # 7,8-dihydroisoquinoline (C5=C6 sp2; C7,C8 sp3; N2 arom → no 2-methyl)
    ("C1=Cc2ccncc2CC1", "7,8-dihydroisoquinoline"),
    ("Cc1nccc2c1CCC=C2", "1-methyl-7,8-dihydroisoquinoline"),
    ("Cc1cc2c(cn1)CCC=C2", "3-methyl-7,8-dihydroisoquinoline"),
    ("Cc1cncc2c1C=CCC2", "4-methyl-7,8-dihydroisoquinoline"),
    ("CC1=CCCc2cnccc21", "5-methyl-7,8-dihydroisoquinoline"),
    ("CC1=Cc2ccncc2CC1", "6-methyl-7,8-dihydroisoquinoline"),
    ("CC1C=Cc2ccncc2C1", "7-methyl-7,8-dihydroisoquinoline"),
    ("CC1CC=Cc2ccncc21", "8-methyl-7,8-dihydroisoquinoline"),
])
def test_phase656(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
