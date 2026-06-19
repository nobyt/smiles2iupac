"""Phase 636: tetrahydronaphthyridines —
5,6,7,8-tetrahydro-1,5/1,6/1,7/1,8-naphthyridine (aromatic pyridine ring, saturated ring with N)
and 1,2,3,4-tetrahydro-1,6/1,7-naphthyridine (aromatic pyridine ring, saturated ring with N at pos 1)
(IUPAC 2013).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("c1cnc2c(c1)NCCC2",  "5,6,7,8-tetrahydro-1,5-naphthyridine"),
    ("c1cnc2c(c1)CNCC2",  "5,6,7,8-tetrahydro-1,6-naphthyridine"),
    ("c1cnc2c(c1)CCNC2",  "5,6,7,8-tetrahydro-1,7-naphthyridine"),
    ("c1cnc2c(c1)CCCN2",  "5,6,7,8-tetrahydro-1,8-naphthyridine"),
    ("c1cc2c(cn1)CCCN2",  "1,2,3,4-tetrahydro-1,6-naphthyridine"),
    ("c1cc2c(cn1)NCCC2",  "1,2,3,4-tetrahydro-1,7-naphthyridine"),
    # 5,6,7,8-tetrahydro-1,5-naphthyridine: N(1,None)-C(2)-C(3)-C(4)-C(4a,junc)-C(8a,junc)-N(5,NH)-C(6)-C(7)-C(8)
    ("Cc1cnc2c(c1)NCCC2",   "3-methyl-5,6,7,8-tetrahydro-1,5-naphthyridine"),
    ("Cc1ccc2c(n1)CCCN2",   "2-methyl-5,6,7,8-tetrahydro-1,5-naphthyridine"),
    ("Cc1ccnc2c1NCCC2",     "4-methyl-5,6,7,8-tetrahydro-1,5-naphthyridine"),
    ("CN1CCCc2ncccc21",     "5-methyl-5,6,7,8-tetrahydro-1,5-naphthyridine"),
    ("CC1CCc2ncccc2N1",     "6-methyl-5,6,7,8-tetrahydro-1,5-naphthyridine"),
    ("CC1CNc2cccnc2C1",     "7-methyl-5,6,7,8-tetrahydro-1,5-naphthyridine"),
    ("CC1CCNc2cccnc21",     "8-methyl-5,6,7,8-tetrahydro-1,5-naphthyridine"),
    # 5,6,7,8-tetrahydro-1,6-naphthyridine: N(1,None)-C(2)-C(3)-C(4)-C(4a,junc)-C(8a,junc)-C(5)-N(6,NH)-C(7)-C(8)
    ("Cc1cnc2c(c1)CNCC2",   "3-methyl-5,6,7,8-tetrahydro-1,6-naphthyridine"),
    ("Cc1ccc2c(n1)CCNC2",   "2-methyl-5,6,7,8-tetrahydro-1,6-naphthyridine"),
    ("Cc1ccnc2c1CNCC2",     "4-methyl-5,6,7,8-tetrahydro-1,6-naphthyridine"),
    ("CC1NCCc2ncccc21",     "5-methyl-5,6,7,8-tetrahydro-1,6-naphthyridine"),
    ("CN1CCc2ncccc2C1",     "6-methyl-5,6,7,8-tetrahydro-1,6-naphthyridine"),
    ("CC1Cc2ncccc2CN1",     "7-methyl-5,6,7,8-tetrahydro-1,6-naphthyridine"),
    ("CC1CNCc2cccnc21",     "8-methyl-5,6,7,8-tetrahydro-1,6-naphthyridine"),
    # 5,6,7,8-tetrahydro-1,7-naphthyridine: N(1,None)-C(2)-C(3)-C(4)-C(4a,junc)-C(8a,junc)-C(5)-C(6)-N(7,NH)-C(8)
    ("Cc1cnc2c(c1)CCNC2",   "3-methyl-5,6,7,8-tetrahydro-1,7-naphthyridine"),
    ("Cc1ccc2c(n1)CNCC2",   "2-methyl-5,6,7,8-tetrahydro-1,7-naphthyridine"),
    ("Cc1ccnc2c1CCNC2",     "4-methyl-5,6,7,8-tetrahydro-1,7-naphthyridine"),
    ("CC1CNCc2ncccc21",     "5-methyl-5,6,7,8-tetrahydro-1,7-naphthyridine"),
    ("CC1Cc2cccnc2CN1",     "6-methyl-5,6,7,8-tetrahydro-1,7-naphthyridine"),
    ("CN1CCc2cccnc2C1",     "7-methyl-5,6,7,8-tetrahydro-1,7-naphthyridine"),
    ("CC1NCCc2cccnc21",     "8-methyl-5,6,7,8-tetrahydro-1,7-naphthyridine"),
    # 5,6,7,8-tetrahydro-1,8-naphthyridine: N(1,None)-C(2)-C(3)-C(4)-C(4a,junc)-C(8a,junc)-C(5)-C(6)-C(7)-N(8,NH)
    ("Cc1cnc2c(c1)CCCN2",   "3-methyl-5,6,7,8-tetrahydro-1,8-naphthyridine"),
    ("Cc1ccc2c(n1)NCCC2",   "2-methyl-5,6,7,8-tetrahydro-1,8-naphthyridine"),
    ("Cc1ccnc2c1CCCN2",     "4-methyl-5,6,7,8-tetrahydro-1,8-naphthyridine"),
    ("CC1CCNc2ncccc21",     "5-methyl-5,6,7,8-tetrahydro-1,8-naphthyridine"),
    ("CC1CNc2ncccc2C1",     "6-methyl-5,6,7,8-tetrahydro-1,8-naphthyridine"),
    ("CC1CCc2cccnc2N1",     "7-methyl-5,6,7,8-tetrahydro-1,8-naphthyridine"),
    ("CN1CCCc2cccnc21",     "8-methyl-5,6,7,8-tetrahydro-1,8-naphthyridine"),
    # 1,2,3,4-tetrahydro-1,6-naphthyridine: N(1,NH)-C(2)-C(3)-C(4)-C(4a,junc)-C(5)-N(6,None)-C(7)-C(8)-C(8a,junc)
    ("CN1CCCc2cnccc21",     "1-methyl-1,2,3,4-tetrahydro-1,6-naphthyridine"),
    ("CC1CCc2cnccc2N1",     "2-methyl-1,2,3,4-tetrahydro-1,6-naphthyridine"),
    ("CC1CNc2ccncc2C1",     "3-methyl-1,2,3,4-tetrahydro-1,6-naphthyridine"),
    ("CC1CCNc2ccncc21",     "4-methyl-1,2,3,4-tetrahydro-1,6-naphthyridine"),
    ("Cc1nccc2c1CCCN2",     "5-methyl-1,2,3,4-tetrahydro-1,6-naphthyridine"),
    ("Cc1cc2c(cn1)CCCN2",   "7-methyl-1,2,3,4-tetrahydro-1,6-naphthyridine"),
    ("Cc1cncc2c1NCCC2",     "8-methyl-1,2,3,4-tetrahydro-1,6-naphthyridine"),
    # 1,2,3,4-tetrahydro-1,7-naphthyridine: N(1,NH)-C(2)-C(3)-C(4)-C(4a,junc)-C(5)-C(6)-N(7,None)-C(8)-C(8a,junc)
    ("CN1CCCc2ccncc21",     "1-methyl-1,2,3,4-tetrahydro-1,7-naphthyridine"),
    ("CC1CCc2ccncc2N1",     "2-methyl-1,2,3,4-tetrahydro-1,7-naphthyridine"),
    ("CC1CNc2cnccc2C1",     "3-methyl-1,2,3,4-tetrahydro-1,7-naphthyridine"),
    ("CC1CCNc2cnccc21",     "4-methyl-1,2,3,4-tetrahydro-1,7-naphthyridine"),
    ("Cc1cncc2c1CCCN2",     "5-methyl-1,2,3,4-tetrahydro-1,7-naphthyridine"),
    ("Cc1nccc2c1NCCC2",     "8-methyl-1,2,3,4-tetrahydro-1,7-naphthyridine"),
    ("Cc1cc2c(cn1)NCCC2",   "6-methyl-1,2,3,4-tetrahydro-1,7-naphthyridine"),
])
def test_phase636_tetrahydronaphthyridines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
