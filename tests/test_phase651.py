"""Phase 651: 1,2-dihydro-1,5-naphthyridine, 3,4-dihydro-1,5-naphthyridine,
3,4-dihydro-2,6-naphthyridine, 1,2-dihydro-2,7-naphthyridine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2-dihydro-1,5-naphthyridine (N5 arom → no 5-methyl)
    ("C1=Cc2ncccc2NC1", "1,2-dihydro-1,5-naphthyridine"),
    ("CN1CC=Cc2ncccc21", "1-methyl-1,2-dihydro-1,5-naphthyridine"),
    ("CC1C=Cc2ncccc2N1", "2-methyl-1,2-dihydro-1,5-naphthyridine"),
    ("CC1=Cc2ncccc2NC1", "3-methyl-1,2-dihydro-1,5-naphthyridine"),
    ("CC1=CCNc2cccnc21", "4-methyl-1,2-dihydro-1,5-naphthyridine"),
    ("Cc1ccc2c(n1)C=CCN2", "6-methyl-1,2-dihydro-1,5-naphthyridine"),
    ("Cc1cnc2c(c1)NCC=C2", "7-methyl-1,2-dihydro-1,5-naphthyridine"),
    ("Cc1ccnc2c1NCC=C2", "8-methyl-1,2-dihydro-1,5-naphthyridine"),
    # 3,4-dihydro-1,5-naphthyridine (N1 sp2 imine → no 1-methyl; N5 arom → no 5-methyl)
    ("C1=Nc2cccnc2CC1", "3,4-dihydro-1,5-naphthyridine"),
    ("CC1=Nc2cccnc2CC1", "2-methyl-3,4-dihydro-1,5-naphthyridine"),
    ("CC1C=Nc2cccnc2C1", "3-methyl-3,4-dihydro-1,5-naphthyridine"),
    ("CC1CC=Nc2cccnc21", "4-methyl-3,4-dihydro-1,5-naphthyridine"),
    ("Cc1ccc2c(n1)CCC=N2", "6-methyl-3,4-dihydro-1,5-naphthyridine"),
    ("Cc1cnc2c(c1)N=CCC2", "7-methyl-3,4-dihydro-1,5-naphthyridine"),
    ("Cc1ccnc2c1N=CCC2", "8-methyl-3,4-dihydro-1,5-naphthyridine"),
    # 3,4-dihydro-2,6-naphthyridine (N2 sp2 imine → no 2-methyl; N6 arom → no 6-methyl)
    ("C1=NCCc2cnccc21", "3,4-dihydro-2,6-naphthyridine"),
    ("CC1=NCCc2cnccc21", "1-methyl-3,4-dihydro-2,6-naphthyridine"),
    ("CC1Cc2cnccc2C=N1", "3-methyl-3,4-dihydro-2,6-naphthyridine"),
    ("CC1CN=Cc2ccncc21", "4-methyl-3,4-dihydro-2,6-naphthyridine"),
    ("Cc1nccc2c1CCN=C2", "5-methyl-3,4-dihydro-2,6-naphthyridine"),
    ("Cc1cc2c(cn1)CCN=C2", "7-methyl-3,4-dihydro-2,6-naphthyridine"),
    ("Cc1cncc2c1C=NCC2", "8-methyl-3,4-dihydro-2,6-naphthyridine"),
    # 1,2-dihydro-2,7-naphthyridine (N7 arom → no 7-methyl)
    ("C1=Cc2ccncc2CN1", "1,2-dihydro-2,7-naphthyridine"),
    ("CC1NC=Cc2ccncc21", "1-methyl-1,2-dihydro-2,7-naphthyridine"),
    ("CN1C=Cc2ccncc2C1", "2-methyl-1,2-dihydro-2,7-naphthyridine"),
    ("CC1=Cc2ccncc2CN1", "3-methyl-1,2-dihydro-2,7-naphthyridine"),
    ("CC1=CNCc2cnccc21", "4-methyl-1,2-dihydro-2,7-naphthyridine"),
    ("Cc1cncc2c1C=CNC2", "5-methyl-1,2-dihydro-2,7-naphthyridine"),
    ("Cc1cc2c(cn1)CNC=C2", "6-methyl-1,2-dihydro-2,7-naphthyridine"),
    ("Cc1nccc2c1CNC=C2", "8-methyl-1,2-dihydro-2,7-naphthyridine"),
])
def test_phase651(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
