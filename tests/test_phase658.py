"""Phase 658: chromane, isochromane, and 1,2,3,4-tetrahydroquinoline methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # chromane (3,4-dihydro-2H-chromene): O(1)-C(2)-C(3)-C(4)-C4a-C(5)-C(6)-C(7)-C(8)-C8a
    # (parent in Phase 133; O1 not methylable)
    ("CC1CCc2ccccc2O1",    "2-methylchromane"),
    ("CC1COc2ccccc2C1",    "3-methylchromane"),
    ("CC1CCOc2ccccc21",    "4-methylchromane"),
    ("Cc1cccc2c1CCCO2",    "5-methylchromane"),
    ("Cc1ccc2c(c1)CCCO2",  "6-methylchromane"),
    ("Cc1ccc2c(c1)OCCC2",  "7-methylchromane"),
    ("Cc1cccc2c1OCCC2",    "8-methylchromane"),
    # isochromane (3,4-dihydro-1H-2-benzopyran): C(1)-O(2)-C(3)-C(4)-C4a-C(5)-C(6)-C(7)-C(8)-C8a
    # (parent + 6-methyl in Phase 406; O2 not methylable)
    ("CC1OCCc2ccccc21",    "1-methylisochromane"),
    ("CC1Cc2ccccc2CO1",    "3-methylisochromane"),
    ("CC1COCc2ccccc21",    "4-methylisochromane"),
    ("Cc1cccc2c1CCOC2",    "5-methylisochromane"),
    ("Cc1ccc2c(c1)CCOC2",  "6-methylisochromane"),
    ("Cc1ccc2c(c1)COCC2",  "7-methylisochromane"),
    ("Cc1cccc2c1COCC2",    "8-methylisochromane"),
    # 1,2,3,4-tetrahydroquinoline: N(1)-C(2)-C(3)-C(4)-C4a-C(5)-C(6)-C(7)-C(8)-C8a
    # (parent in Phase 133; 1-methyl in Phase 404)
    ("CC1CCc2ccccc2N1",    "2-methyl-1,2,3,4-tetrahydroquinoline"),
    ("CC1CNc2ccccc2C1",    "3-methyl-1,2,3,4-tetrahydroquinoline"),
    ("CC1CCNc2ccccc21",    "4-methyl-1,2,3,4-tetrahydroquinoline"),
    ("Cc1cccc2c1CCCN2",    "5-methyl-1,2,3,4-tetrahydroquinoline"),
    ("Cc1ccc2c(c1)CCCN2",  "6-methyl-1,2,3,4-tetrahydroquinoline"),
    ("Cc1ccc2c(c1)NCCC2",  "7-methyl-1,2,3,4-tetrahydroquinoline"),
    ("Cc1cccc2c1NCCC2",    "8-methyl-1,2,3,4-tetrahydroquinoline"),
])
def test_phase658(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
