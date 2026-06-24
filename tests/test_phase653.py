"""Phase 653: 1,4-dihydroquinoline, 1,4-dihydroisoquinoline, 1,2,3,4-tetrahydrocinnoline."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,4-dihydroquinoline (N1H sp2 enamine; C2=C3; C4 sp3; benzo)
    ("C1=CNc2ccccc2C1", "1,4-dihydroquinoline"),
    ("CN1C=CCc2ccccc21", "1-methyl-1,4-dihydroquinoline"),
    ("CC1=CNc2ccccc2C1", "2-methyl-1,4-dihydroquinoline"),
    ("CC1=CCc2ccccc2N1", "3-methyl-1,4-dihydroquinoline"),
    ("CC1C=CNc2ccccc21", "4-methyl-1,4-dihydroquinoline"),
    ("Cc1cccc2c1CC=CN2", "5-methyl-1,4-dihydroquinoline"),
    ("Cc1ccc2c(c1)CC=CN2", "6-methyl-1,4-dihydroquinoline"),
    ("Cc1ccc2c(c1)NC=CC2", "7-methyl-1,4-dihydroquinoline"),
    ("Cc1cccc2c1NC=CC2", "8-methyl-1,4-dihydroquinoline"),
    # 1,4-dihydroisoquinoline (N2 imine no H; C1H2, C4H2 sp3; C3H sp2; benzo)
    ("C1=NCc2ccccc2C1", "1,4-dihydroisoquinoline"),
    ("CC1N=CCc2ccccc21", "1-methyl-1,4-dihydroisoquinoline"),
    ("CC1=NCc2ccccc2C1", "3-methyl-1,4-dihydroisoquinoline"),
    ("CC1C=NCc2ccccc21", "4-methyl-1,4-dihydroisoquinoline"),
    ("Cc1cccc2c1CC=NC2", "5-methyl-1,4-dihydroisoquinoline"),
    ("Cc1ccc2c(c1)CC=NC2", "6-methyl-1,4-dihydroisoquinoline"),
    ("Cc1ccc2c(c1)CN=CC2", "7-methyl-1,4-dihydroisoquinoline"),
    ("Cc1cccc2c1CN=CC2", "8-methyl-1,4-dihydroisoquinoline"),
    # 1,2,3,4-tetrahydrocinnoline (N1H sp2; N2H sp3; C3H2, C4H2 sp3; benzo)
    ("c1ccc2c(c1)CCNN2", "1,2,3,4-tetrahydrocinnoline"),
    ("CN1NCCc2ccccc21", "1-methyl-1,2,3,4-tetrahydrocinnoline"),
    ("CN1CCc2ccccc2N1", "2-methyl-1,2,3,4-tetrahydrocinnoline"),
    ("CC1Cc2ccccc2NN1", "3-methyl-1,2,3,4-tetrahydrocinnoline"),
    ("CC1CNNc2ccccc21", "4-methyl-1,2,3,4-tetrahydrocinnoline"),
    ("Cc1cccc2c1CCNN2", "5-methyl-1,2,3,4-tetrahydrocinnoline"),
    ("Cc1ccc2c(c1)CCNN2", "6-methyl-1,2,3,4-tetrahydrocinnoline"),
    ("Cc1ccc2c(c1)NNCC2", "7-methyl-1,2,3,4-tetrahydrocinnoline"),
    ("Cc1cccc2c1NNCC2", "8-methyl-1,2,3,4-tetrahydrocinnoline"),
])
def test_phase653(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
