"""Phase 685: chromane and isochromane methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # chromane: O(1)-C(2)-C(3)-C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent retained name; O1 not methylable; C2-C4 sp3 and C5-C8 aromatic methylable)
    ("c1ccc2c(c1)CCCO2",   "chromane"),
    ("CC1CCc2ccccc2O1",    "2-methylchromane"),
    ("CC1COc2ccccc2C1",    "3-methylchromane"),
    ("CC1CCOc2ccccc21",    "4-methylchromane"),
    ("Cc1cccc2c1CCCO2",    "5-methylchromane"),
    ("Cc1ccc2c(c1)CCCO2",  "6-methylchromane"),
    ("Cc1ccc2c(c1)OCCC2",  "7-methylchromane"),
    ("Cc1cccc2c1OCCC2",    "8-methylchromane"),
    # isochromane: C(1)-O(2)-C(3)-C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent retained name; O2 not methylable; C1, C3-C4 sp3 and C5-C8 aromatic methylable)
    ("c1ccc2c(c1)CCOC2",   "isochromane"),
    ("CC1OCCc2ccccc21",    "1-methylisochromane"),
    ("CC1Cc2ccccc2CO1",    "3-methylisochromane"),
    ("CC1COCc2ccccc21",    "4-methylisochromane"),
    ("Cc1cccc2c1CCOC2",    "5-methylisochromane"),
    ("Cc1ccc2c(c1)CCOC2",  "6-methylisochromane"),
    ("Cc1ccc2c(c1)COCC2",  "7-methylisochromane"),
    ("Cc1cccc2c1COCC2",    "8-methylisochromane"),
])
def test_phase685(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
