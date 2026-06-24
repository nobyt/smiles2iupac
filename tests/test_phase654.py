"""Phase 654: fix 3,4-dihydro-1,8-/2,7-naphthyridine (C2-sym rename); add 7,8-dihydro-1,7-naphthyridine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 7,8-dihydro-1,7-naphthyridine (N7H sp2 enamine; C8 sp3; C5=C6; N1 arom)
    ("C1=Cc2cccnc2CN1", "7,8-dihydro-1,7-naphthyridine"),
    ("Cc1ccc2c(n1)CNC=C2", "2-methyl-7,8-dihydro-1,7-naphthyridine"),
    ("Cc1cnc2c(c1)C=CNC2", "3-methyl-7,8-dihydro-1,7-naphthyridine"),
    ("Cc1ccnc2c1C=CNC2", "4-methyl-7,8-dihydro-1,7-naphthyridine"),
    ("CC1=CNCc2ncccc21", "5-methyl-7,8-dihydro-1,7-naphthyridine"),
    ("CC1=Cc2cccnc2CN1", "6-methyl-7,8-dihydro-1,7-naphthyridine"),
    ("CN1C=Cc2cccnc2C1", "7-methyl-7,8-dihydro-1,7-naphthyridine"),
    ("CC1NC=Cc2cccnc21", "8-methyl-7,8-dihydro-1,7-naphthyridine"),
])
def test_phase654(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
