"""Phase 652: 1,2-dihydrocinnoline and 3,4-dihydrocinnoline + methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2-dihydrocinnoline (all locants 1-8 have H)
    ("C1=Cc2ccccc2NN1", "1,2-dihydrocinnoline"),
    ("CN1NC=Cc2ccccc21", "1-methyl-1,2-dihydrocinnoline"),
    ("CN1C=Cc2ccccc2N1", "2-methyl-1,2-dihydrocinnoline"),
    ("CC1=Cc2ccccc2NN1", "3-methyl-1,2-dihydrocinnoline"),
    ("CC1=CNNc2ccccc21", "4-methyl-1,2-dihydrocinnoline"),
    ("Cc1cccc2c1C=CNN2", "5-methyl-1,2-dihydrocinnoline"),
    ("Cc1ccc2c(c1)C=CNN2", "6-methyl-1,2-dihydrocinnoline"),
    ("Cc1ccc2c(c1)NNC=C2", "7-methyl-1,2-dihydrocinnoline"),
    ("Cc1cccc2c1NNC=C2", "8-methyl-1,2-dihydrocinnoline"),
    # 3,4-dihydrocinnoline (N1=N2 azo sp2 → no 1- or 2-methyl)
    ("c1ccc2c(c1)CCN=N2", "3,4-dihydrocinnoline"),
    ("CC1Cc2ccccc2N=N1", "3-methyl-3,4-dihydrocinnoline"),
    ("CC1CN=Nc2ccccc21", "4-methyl-3,4-dihydrocinnoline"),
    ("Cc1cccc2c1CCN=N2", "5-methyl-3,4-dihydrocinnoline"),
    ("Cc1ccc2c(c1)CCN=N2", "6-methyl-3,4-dihydrocinnoline"),
    ("Cc1ccc2c(c1)N=NCC2", "7-methyl-3,4-dihydrocinnoline"),
    ("Cc1cccc2c1N=NCC2", "8-methyl-3,4-dihydrocinnoline"),
])
def test_phase652(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
