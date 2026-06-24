"""Phase 655: 7,8-dihydro-1,6-naphthyridine (N6 sp2 imine; C7,C8 sp3; N1 arom)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 7,8-dihydro-1,6-naphthyridine (C5=N6 imine; C7,C8 sp3; N1 arom → no 1-methyl; N6 imine → no 6-methyl)
    ("C1=NCCc2ncccc21", "7,8-dihydro-1,6-naphthyridine"),
    ("Cc1ccc2c(n1)CCN=C2", "2-methyl-7,8-dihydro-1,6-naphthyridine"),
    ("Cc1cnc2c(c1)C=NCC2", "3-methyl-7,8-dihydro-1,6-naphthyridine"),
    ("Cc1ccnc2c1C=NCC2", "4-methyl-7,8-dihydro-1,6-naphthyridine"),
    ("CC1=NCCc2ncccc21", "5-methyl-7,8-dihydro-1,6-naphthyridine"),
    ("CC1Cc2ncccc2C=N1", "7-methyl-7,8-dihydro-1,6-naphthyridine"),
    ("CC1CN=Cc2cccnc21", "8-methyl-7,8-dihydro-1,6-naphthyridine"),
])
def test_phase655(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
