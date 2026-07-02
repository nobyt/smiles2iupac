"""Phase 685: chromane and isochromane methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # chromane: O(1)-C(2)-C(3)-C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent retained name; O1 not methylable; C2-C4 sp3 and C5-C8 aromatic methylable)
    ("c1ccc2c(c1)CCCO2",   "3,4-dihydro-2H-chromene"),
    ("CC1CCc2ccccc2O1",    "2-methyl-3,4-dihydro-2H-chromene"),
    ("CC1COc2ccccc2C1",    "3-methyl-3,4-dihydro-2H-chromene"),
    ("CC1CCOc2ccccc21",    "4-methyl-3,4-dihydro-2H-chromene"),
    ("Cc1cccc2c1CCCO2",    "5-methyl-3,4-dihydro-2H-chromene"),
    ("Cc1ccc2c(c1)CCCO2",  "6-methyl-3,4-dihydro-2H-chromene"),
    ("Cc1ccc2c(c1)OCCC2",  "7-methyl-3,4-dihydro-2H-chromene"),
    ("Cc1cccc2c1OCCC2",    "8-methyl-3,4-dihydro-2H-chromene"),
    # isochromane: C(1)-O(2)-C(3)-C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent retained name; O2 not methylable; C1, C3-C4 sp3 and C5-C8 aromatic methylable)
    ("c1ccc2c(c1)CCOC2",   "3,4-dihydro-1H-isochromene"),
    ("CC1OCCc2ccccc21",    "1-methyl-3,4-dihydro-1H-isochromene"),
    ("CC1Cc2ccccc2CO1",    "3-methyl-3,4-dihydro-1H-isochromene"),
    ("CC1COCc2ccccc21",    "4-methyl-3,4-dihydro-1H-isochromene"),
    ("Cc1cccc2c1CCOC2",    "5-methyl-3,4-dihydro-1H-isochromene"),
    ("Cc1ccc2c(c1)CCOC2",  "6-methyl-3,4-dihydro-1H-isochromene"),
    ("Cc1ccc2c(c1)COCC2",  "7-methyl-3,4-dihydro-1H-isochromene"),
    ("Cc1cccc2c1COCC2",    "8-methyl-3,4-dihydro-1H-isochromene"),
])
def test_phase685(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
