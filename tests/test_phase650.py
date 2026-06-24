"""Phase 650: 3,4-dihydro-1,6-naphthyridine, 5,6-dihydro-1,7-naphthyridine,
1,2-dihydro-1,7-naphthyridine, 3,4-dihydro-1,7-naphthyridine,
3,4-dihydro-1,8-naphthyridine, 3,4-dihydro-2,7-naphthyridine,
1,2-dihydro-2,6-naphthyridine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 3,4-dihydro-1,6-naphthyridine (OPSIN-verified; N1 sp2 imine → no 1-methyl; N6 arom → no 6-methyl)
    ("C1=Nc2ccncc2CC1", "3,4-dihydro-1,6-naphthyridine"),
    ("CC1=Nc2ccncc2CC1", "2-methyl-3,4-dihydro-1,6-naphthyridine"),
    ("CC1C=Nc2ccncc2C1", "3-methyl-3,4-dihydro-1,6-naphthyridine"),
    ("CC1CC=Nc2ccncc21", "4-methyl-3,4-dihydro-1,6-naphthyridine"),
    ("Cc1nccc2c1CCC=N2", "5-methyl-3,4-dihydro-1,6-naphthyridine"),
    ("Cc1cc2c(cn1)CCC=N2", "7-methyl-3,4-dihydro-1,6-naphthyridine"),
    ("Cc1cncc2c1N=CCC2", "8-methyl-3,4-dihydro-1,6-naphthyridine"),
    # 5,6-dihydro-1,7-naphthyridine (OPSIN-verified; N1 sp2 imine → no 1-methyl; N7 arom → no 7-methyl)
    ("C1=NCCc2cccnc21", "5,6-dihydro-1,7-naphthyridine"),
    ("Cc1ccc2c(n1)C=NCC2", "2-methyl-5,6-dihydro-1,7-naphthyridine"),
    ("Cc1cnc2c(c1)CCN=C2", "3-methyl-5,6-dihydro-1,7-naphthyridine"),
    ("Cc1ccnc2c1CCN=C2", "4-methyl-5,6-dihydro-1,7-naphthyridine"),
    ("CC1CN=Cc2ncccc21", "5-methyl-5,6-dihydro-1,7-naphthyridine"),
    ("CC1Cc2cccnc2C=N1", "6-methyl-5,6-dihydro-1,7-naphthyridine"),
    ("CC1=NCCc2cccnc21", "8-methyl-5,6-dihydro-1,7-naphthyridine"),
    # 1,2-dihydro-1,7-naphthyridine (OPSIN-verified; N7 arom → no 7-methyl)
    ("C1=Cc2ccncc2NC1", "1,2-dihydro-1,7-naphthyridine"),
    ("CN1CC=Cc2ccncc21", "1-methyl-1,2-dihydro-1,7-naphthyridine"),
    ("CC1C=Cc2ccncc2N1", "2-methyl-1,2-dihydro-1,7-naphthyridine"),
    ("CC1=Cc2ccncc2NC1", "3-methyl-1,2-dihydro-1,7-naphthyridine"),
    ("CC1=CCNc2cnccc21", "4-methyl-1,2-dihydro-1,7-naphthyridine"),
    ("Cc1cncc2c1C=CCN2", "5-methyl-1,2-dihydro-1,7-naphthyridine"),
    ("Cc1cc2c(cn1)NCC=C2", "6-methyl-1,2-dihydro-1,7-naphthyridine"),
    ("Cc1nccc2c1NCC=C2", "8-methyl-1,2-dihydro-1,7-naphthyridine"),
    # 3,4-dihydro-1,7-naphthyridine (OPSIN-verified; N1 sp2 imine → no 1-methyl; N7 arom → no 7-methyl)
    ("C1=Nc2cnccc2CC1", "3,4-dihydro-1,7-naphthyridine"),
    ("CC1=Nc2cnccc2CC1", "2-methyl-3,4-dihydro-1,7-naphthyridine"),
    ("CC1C=Nc2cnccc2C1", "3-methyl-3,4-dihydro-1,7-naphthyridine"),
    ("CC1CC=Nc2cnccc21", "4-methyl-3,4-dihydro-1,7-naphthyridine"),
    ("Cc1cncc2c1CCC=N2", "5-methyl-3,4-dihydro-1,7-naphthyridine"),
    ("Cc1cc2c(cn1)N=CCC2", "6-methyl-3,4-dihydro-1,7-naphthyridine"),
    ("Cc1nccc2c1N=CCC2", "8-methyl-3,4-dihydro-1,7-naphthyridine"),
    # 3,4-dihydro-1,8-naphthyridine (C2-sym rename; N1 arom → no 1-methyl; N8 sp2 imine → no 8-methyl)
    ("C1=Nc2ncccc2CC1", "3,4-dihydro-1,8-naphthyridine"),
    ("CC1=Nc2ncccc2CC1", "2-methyl-3,4-dihydro-1,8-naphthyridine"),
    ("CC1C=Nc2ncccc2C1", "3-methyl-3,4-dihydro-1,8-naphthyridine"),
    ("CC1CC=Nc2ncccc21", "4-methyl-3,4-dihydro-1,8-naphthyridine"),
    ("Cc1ccnc2c1CCC=N2", "5-methyl-3,4-dihydro-1,8-naphthyridine"),
    ("Cc1cnc2c(c1)CCC=N2", "6-methyl-3,4-dihydro-1,8-naphthyridine"),
    ("Cc1ccc2c(n1)N=CCC2", "7-methyl-3,4-dihydro-1,8-naphthyridine"),
    # 3,4-dihydro-2,7-naphthyridine (C2-sym rename; N2 arom → no 2-methyl; N7 sp2 imine → no 7-methyl)
    ("C1=NCCc2ccncc21", "3,4-dihydro-2,7-naphthyridine"),
    ("CC1=NCCc2ccncc21", "1-methyl-3,4-dihydro-2,7-naphthyridine"),
    ("CC1Cc2ccncc2C=N1", "3-methyl-3,4-dihydro-2,7-naphthyridine"),
    ("CC1CN=Cc2cnccc21", "4-methyl-3,4-dihydro-2,7-naphthyridine"),
    ("Cc1cncc2c1CCN=C2", "5-methyl-3,4-dihydro-2,7-naphthyridine"),
    ("Cc1cc2c(cn1)C=NCC2", "6-methyl-3,4-dihydro-2,7-naphthyridine"),
    ("Cc1nccc2c1C=NCC2", "8-methyl-3,4-dihydro-2,7-naphthyridine"),
    # 1,2-dihydro-2,6-naphthyridine (OPSIN-verified; N6 arom → no 6-methyl)
    ("C1=Cc2cnccc2CN1", "1,2-dihydro-2,6-naphthyridine"),
    ("CC1NC=Cc2cnccc21", "1-methyl-1,2-dihydro-2,6-naphthyridine"),
    ("CN1C=Cc2cnccc2C1", "2-methyl-1,2-dihydro-2,6-naphthyridine"),
    ("CC1=Cc2cnccc2CN1", "3-methyl-1,2-dihydro-2,6-naphthyridine"),
    ("CC1=CNCc2ccncc21", "4-methyl-1,2-dihydro-2,6-naphthyridine"),
    ("Cc1nccc2c1C=CNC2", "5-methyl-1,2-dihydro-2,6-naphthyridine"),
    ("Cc1cc2c(cn1)C=CNC2", "7-methyl-1,2-dihydro-2,6-naphthyridine"),
    ("Cc1cncc2c1CNC=C2", "8-methyl-1,2-dihydro-2,6-naphthyridine"),
])
def test_phase650(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
