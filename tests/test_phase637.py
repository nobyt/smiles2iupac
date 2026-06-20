"""Phase 637: 1,2,3,4-tetrahydro-2,6-naphthyridine and 1,2,3,4-tetrahydro-2,7-naphthyridine —
symmetric naphthyridines (N at 2 and 6, or 2 and 7); 1,2,3,4-tetrahydro prefix is preferred
since {1,2,3,4} < {5,6,7,8} (IUPAC 2013 lowest-locant rule).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("c1cc2c(cn1)CCNC2",  "1,2,3,4-tetrahydro-2,6-naphthyridine"),
    ("c1cc2c(cn1)CNCC2",  "1,2,3,4-tetrahydro-2,7-naphthyridine"),
    # 1,2,3,4-tetrahydro-2,6-naphthyridine: C(1)-N(2,H)-C(3)-C(4)-C(4a,junc)-C(5)-N(6,None)-C(7)-C(8)-C(8a,junc)
    ("CC1NCCc2cnccc21",     "1-methyl-1,2,3,4-tetrahydro-2,6-naphthyridine"),
    ("CN1CCc2cnccc2C1",     "2-methyl-1,2,3,4-tetrahydro-2,6-naphthyridine"),
    ("CC1Cc2cnccc2CN1",     "3-methyl-1,2,3,4-tetrahydro-2,6-naphthyridine"),
    ("CC1CNCc2ccncc21",     "4-methyl-1,2,3,4-tetrahydro-2,6-naphthyridine"),
    ("Cc1nccc2c1CCNC2",     "5-methyl-1,2,3,4-tetrahydro-2,6-naphthyridine"),
    ("Cc1cc2c(cn1)CCNC2",   "7-methyl-1,2,3,4-tetrahydro-2,6-naphthyridine"),
    ("Cc1cncc2c1CNCC2",     "8-methyl-1,2,3,4-tetrahydro-2,6-naphthyridine"),
    # 1,2,3,4-tetrahydro-2,7-naphthyridine: C(1)-N(2,H)-C(3)-C(4)-C(4a,junc)-C(5)-C(6)-N(7,None)-C(8)-C(8a,junc)
    ("CC1NCCc2ccncc21",     "1-methyl-1,2,3,4-tetrahydro-2,7-naphthyridine"),
    ("CN1CCc2ccncc2C1",     "2-methyl-1,2,3,4-tetrahydro-2,7-naphthyridine"),
    ("CC1Cc2ccncc2CN1",     "3-methyl-1,2,3,4-tetrahydro-2,7-naphthyridine"),
    ("CC1CNCc2cnccc21",     "4-methyl-1,2,3,4-tetrahydro-2,7-naphthyridine"),
    ("Cc1cncc2c1CCNC2",     "5-methyl-1,2,3,4-tetrahydro-2,7-naphthyridine"),
    ("Cc1cc2c(cn1)CNCC2",   "6-methyl-1,2,3,4-tetrahydro-2,7-naphthyridine"),
    ("Cc1nccc2c1CNCC2",     "8-methyl-1,2,3,4-tetrahydro-2,7-naphthyridine"),
])
def test_phase637_tetrahydronaphthyridines_2x(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
