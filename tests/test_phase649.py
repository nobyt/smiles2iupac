"""Phase 649: 1,2-dihydroisoquinoline, 1,2-dihydronaphthyridine,
5,6-dihydro-1,6-naphthyridine, 1,2-dihydro-1,6-naphthyridine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2-dihydroisoquinoline parent + methyl derivatives (OPSIN-verified; locants 1-8)
    ("C1=Cc2ccccc2CN1", "1,2-dihydroisoquinoline"),
    ("CC1NC=Cc2ccccc21", "1-methyl-1,2-dihydroisoquinoline"),
    ("CN1C=Cc2ccccc2C1", "2-methyl-1,2-dihydroisoquinoline"),
    ("CC1=Cc2ccccc2CN1", "3-methyl-1,2-dihydroisoquinoline"),
    ("CC1=CNCc2ccccc21", "4-methyl-1,2-dihydroisoquinoline"),
    ("Cc1cccc2c1C=CNC2", "5-methyl-1,2-dihydroisoquinoline"),
    ("Cc1ccc2c(c1)C=CNC2", "6-methyl-1,2-dihydroisoquinoline"),
    ("Cc1ccc2c(c1)CNC=C2", "7-methyl-1,2-dihydroisoquinoline"),
    ("Cc1cccc2c1CNC=C2", "8-methyl-1,2-dihydroisoquinoline"),
    # 1,2-dihydronaphthyridine parent + methyl derivatives (OPSIN-verified; N8 arom → no 8-methyl)
    ("C1=Cc2cccnc2NC1", "1,2-dihydronaphthyridine"),
    ("CN1CC=Cc2cccnc21", "1-methyl-1,2-dihydronaphthyridine"),
    ("CC1C=Cc2cccnc2N1", "2-methyl-1,2-dihydronaphthyridine"),
    ("CC1=Cc2cccnc2NC1", "3-methyl-1,2-dihydronaphthyridine"),
    ("CC1=CCNc2ncccc21", "4-methyl-1,2-dihydronaphthyridine"),
    ("Cc1ccnc2c1C=CCN2", "5-methyl-1,2-dihydronaphthyridine"),
    ("Cc1cnc2c(c1)C=CCN2", "6-methyl-1,2-dihydronaphthyridine"),
    ("Cc1ccc2c(n1)NCC=C2", "7-methyl-1,2-dihydronaphthyridine"),
    # 5,6-dihydro-1,6-naphthyridine parent + methyl derivatives (OPSIN-verified; N1 arom → no 1-methyl)
    ("C1=Cc2ncccc2CN1", "5,6-dihydro-1,6-naphthyridine"),
    ("Cc1ccc2c(n1)C=CNC2", "2-methyl-5,6-dihydro-1,6-naphthyridine"),
    ("Cc1cnc2c(c1)CNC=C2", "3-methyl-5,6-dihydro-1,6-naphthyridine"),
    ("Cc1ccnc2c1CNC=C2", "4-methyl-5,6-dihydro-1,6-naphthyridine"),
    ("CC1NC=Cc2ncccc21", "5-methyl-5,6-dihydro-1,6-naphthyridine"),
    ("CN1C=Cc2ncccc2C1", "6-methyl-5,6-dihydro-1,6-naphthyridine"),
    ("CC1=Cc2ncccc2CN1", "7-methyl-5,6-dihydro-1,6-naphthyridine"),
    ("CC1=CNCc2cccnc21", "8-methyl-5,6-dihydro-1,6-naphthyridine"),
    # 1,2-dihydro-1,6-naphthyridine parent + methyl derivatives (OPSIN-verified; N6 arom → no 6-methyl)
    ("C1=Cc2cnccc2NC1", "1,2-dihydro-1,6-naphthyridine"),
    ("CN1CC=Cc2cnccc21", "1-methyl-1,2-dihydro-1,6-naphthyridine"),
    ("CC1C=Cc2cnccc2N1", "2-methyl-1,2-dihydro-1,6-naphthyridine"),
    ("CC1=Cc2cnccc2NC1", "3-methyl-1,2-dihydro-1,6-naphthyridine"),
    ("CC1=CCNc2ccncc21", "4-methyl-1,2-dihydro-1,6-naphthyridine"),
    ("Cc1nccc2c1C=CCN2", "5-methyl-1,2-dihydro-1,6-naphthyridine"),
    ("Cc1cc2c(cn1)C=CCN2", "7-methyl-1,2-dihydro-1,6-naphthyridine"),
    ("Cc1cncc2c1NCC=C2", "8-methyl-1,2-dihydro-1,6-naphthyridine"),
])
def test_phase649(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
