"""Phase 640: 9,10-dihydrophenanthridine and 5,6-dihydrophenanthridine —
partially saturated tricyclic phenanthridine derivatives (IUPAC 2013).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("C1=Cc2cnc3ccccc3c2CC1",    "9,10-dihydrophenanthridine"),
    ("c1ccc2c(c1)CNc1ccccc1-2",  "5,6-dihydrophenanthridine"),
    # 9,10-dihydrophenanthridine: N5(arom)->None; C9,C10 sp3; C7=C8 double bond retained
    # locant 5 (aromatic N) is not C-substitutable; 9 positions tested
    ("Cc1cccc2ncc3c(c12)CCC=C3",  "1-methyl-9,10-dihydrophenanthridine"),
    ("Cc1ccc2ncc3c(c2c1)CCC=C3",  "2-methyl-9,10-dihydrophenanthridine"),
    ("Cc1ccc2c3c(cnc2c1)C=CCC3",  "3-methyl-9,10-dihydrophenanthridine"),
    ("Cc1cccc2c3c(cnc12)C=CCC3",  "4-methyl-9,10-dihydrophenanthridine"),
    ("Cc1nc2ccccc2c2c1C=CCC2",    "6-methyl-9,10-dihydrophenanthridine"),
    ("CC1=CCCc2c1cnc1ccccc21",    "7-methyl-9,10-dihydrophenanthridine"),
    ("CC1=Cc2cnc3ccccc3c2CC1",    "8-methyl-9,10-dihydrophenanthridine"),
    ("CC1C=Cc2cnc3ccccc3c2C1",    "9-methyl-9,10-dihydrophenanthridine"),
    ("CC1CC=Cc2cnc3ccccc3c21",    "10-methyl-9,10-dihydrophenanthridine"),
    # 5,6-dihydrophenanthridine: N5(sp3,NH); C6(sp3,CH2); all 10 positions tested
    ("Cc1cccc2c1-c1ccccc1CN2",    "1-methyl-5,6-dihydrophenanthridine"),
    ("Cc1ccc2c(c1)-c1ccccc1CN2",  "2-methyl-5,6-dihydrophenanthridine"),
    ("Cc1ccc2c(c1)NCc1ccccc1-2",  "3-methyl-5,6-dihydrophenanthridine"),
    ("Cc1cccc2c1NCc1ccccc1-2",    "4-methyl-5,6-dihydrophenanthridine"),
    ("CN1Cc2ccccc2-c2ccccc21",    "5-methyl-5,6-dihydrophenanthridine"),
    ("CC1Nc2ccccc2-c2ccccc21",    "6-methyl-5,6-dihydrophenanthridine"),
    ("Cc1cccc2c1CNc1ccccc1-2",    "7-methyl-5,6-dihydrophenanthridine"),
    ("Cc1ccc2c(c1)CNc1ccccc1-2",  "8-methyl-5,6-dihydrophenanthridine"),
    ("Cc1ccc2c(c1)-c1ccccc1NC2",  "9-methyl-5,6-dihydrophenanthridine"),
    ("Cc1cccc2c1-c1ccccc1NC2",    "10-methyl-5,6-dihydrophenanthridine"),
])
def test_phase640_dihydrophenanthridines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
