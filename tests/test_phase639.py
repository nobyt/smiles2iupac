"""Phase 639: 1,2,3,4-tetrahydrobenzo[g]quinoline, 1,2,3,4-tetrahydrobenzo[f]quinoline,
and 1,2,3,4-tetrahydrobenzo[h]quinoline — tricyclic benzo-fused tetrahydroquinolines
(IUPAC 2013).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("c1ccc2cc3c(cc2c1)CCCN3",  "1,2,3,4-tetrahydrobenzo[g]quinoline"),
    ("c1ccc2c3c(ccc2c1)NCCC3",  "1,2,3,4-tetrahydrobenzo[f]quinoline"),
    ("c1ccc2c3c(ccc2c1)CCCN3",  "1,2,3,4-tetrahydrobenzo[h]quinoline"),
    # 1,2,3,4-tetrahydrobenzo[g]quinoline: N(1)-C2-C3-C4-C(4a)-C5-C(5a)-C6-C7-C8-C9-C(9a)-C10-C(10a)
    ("CN1CCCc2cc3ccccc3cc21",    "1-methyl-1,2,3,4-tetrahydrobenzo[g]quinoline"),
    ("CC1CCc2cc3ccccc3cc2N1",    "2-methyl-1,2,3,4-tetrahydrobenzo[g]quinoline"),
    ("CC1CNc2cc3ccccc3cc2C1",    "3-methyl-1,2,3,4-tetrahydrobenzo[g]quinoline"),
    ("CC1CCNc2cc3ccccc3cc21",    "4-methyl-1,2,3,4-tetrahydrobenzo[g]quinoline"),
    ("Cc1c2c(cc3ccccc13)NCCC2",  "5-methyl-1,2,3,4-tetrahydrobenzo[g]quinoline"),
    ("Cc1cccc2cc3c(cc12)CCCN3",  "6-methyl-1,2,3,4-tetrahydrobenzo[g]quinoline"),
    ("Cc1ccc2cc3c(cc2c1)CCCN3",  "7-methyl-1,2,3,4-tetrahydrobenzo[g]quinoline"),
    ("Cc1ccc2cc3c(cc2c1)NCCC3",  "8-methyl-1,2,3,4-tetrahydrobenzo[g]quinoline"),
    ("Cc1cccc2cc3c(cc12)NCCC3",  "9-methyl-1,2,3,4-tetrahydrobenzo[g]quinoline"),
    ("Cc1c2c(cc3ccccc13)CCCN2",  "10-methyl-1,2,3,4-tetrahydrobenzo[g]quinoline"),
    # 1,2,3,4-tetrahydrobenzo[f]quinoline: C1-C2-C3-N(4)-C(4a)-C5-C6-C7-C8-C(8a)-C9-C10-C(10a)-C(10b)
    ("CC1CCNc2ccc3ccccc3c21",    "1-methyl-1,2,3,4-tetrahydrobenzo[f]quinoline"),
    ("CC1CNc2ccc3ccccc3c2C1",    "2-methyl-1,2,3,4-tetrahydrobenzo[f]quinoline"),
    ("CC1CCc2c(ccc3ccccc23)N1",  "3-methyl-1,2,3,4-tetrahydrobenzo[f]quinoline"),
    ("CN1CCCc2c1ccc1ccccc21",    "4-methyl-1,2,3,4-tetrahydrobenzo[f]quinoline"),
    ("Cc1cc2ccccc2c2c1NCCC2",    "5-methyl-1,2,3,4-tetrahydrobenzo[f]quinoline"),
    ("Cc1cc2c(c3ccccc13)CCCN2",  "6-methyl-1,2,3,4-tetrahydrobenzo[f]quinoline"),
    ("Cc1cccc2c3c(ccc12)NCCC3",  "7-methyl-1,2,3,4-tetrahydrobenzo[f]quinoline"),
    ("Cc1ccc2c3c(ccc2c1)NCCC3",  "8-methyl-1,2,3,4-tetrahydrobenzo[f]quinoline"),
    ("Cc1ccc2ccc3c(c2c1)CCCN3",  "9-methyl-1,2,3,4-tetrahydrobenzo[f]quinoline"),
    ("Cc1cccc2ccc3c(c12)CCCN3",  "10-methyl-1,2,3,4-tetrahydrobenzo[f]quinoline"),
    # 1,2,3,4-tetrahydrobenzo[h]quinoline: N(1)-C2-C3-C4-C(4a)-C5-C6-C7-C8-C(8a)-C9-C10-C(10a)-C(10b)
    ("CN1CCCc2ccc3ccccc3c21",    "1-methyl-1,2,3,4-tetrahydrobenzo[h]quinoline"),
    ("CC1CCc2ccc3ccccc3c2N1",    "2-methyl-1,2,3,4-tetrahydrobenzo[h]quinoline"),
    ("CC1CNc2c(ccc3ccccc23)C1",  "3-methyl-1,2,3,4-tetrahydrobenzo[h]quinoline"),
    ("CC1CCNc2c1ccc1ccccc21",    "4-methyl-1,2,3,4-tetrahydrobenzo[h]quinoline"),
    ("Cc1cc2ccccc2c2c1CCCN2",    "5-methyl-1,2,3,4-tetrahydrobenzo[h]quinoline"),
    ("Cc1cc2c(c3ccccc13)NCCC2",  "6-methyl-1,2,3,4-tetrahydrobenzo[h]quinoline"),
    ("Cc1cccc2c3c(ccc12)CCCN3",  "7-methyl-1,2,3,4-tetrahydrobenzo[h]quinoline"),
    ("Cc1ccc2c3c(ccc2c1)CCCN3",  "8-methyl-1,2,3,4-tetrahydrobenzo[h]quinoline"),
    ("Cc1ccc2ccc3c(c2c1)NCCC3",  "9-methyl-1,2,3,4-tetrahydrobenzo[h]quinoline"),
    ("Cc1cccc2ccc3c(c12)NCCC3",  "10-methyl-1,2,3,4-tetrahydrobenzo[h]quinoline"),
])
def test_phase639_tetrahydrobenzo_quinolines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
