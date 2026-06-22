"""Phase 641: 1,2,3,4-tetrahydrobenzo[f]isoquinoline, 1,2,3,4-tetrahydrobenzo[g]isoquinoline,
and 1,2,3,4-tetrahydrobenzo[h]isoquinoline — tricyclic benzo-fused tetrahydroisoquinolines
(IUPAC 2013).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("c1ccc2c3c(ccc2c1)CNCC3",  "1,2,3,4-tetrahydrobenzo[f]isoquinoline"),
    ("c1ccc2cc3c(cc2c1)CCNC3",  "1,2,3,4-tetrahydrobenzo[g]isoquinoline"),
    ("c1ccc2c3c(ccc2c1)CCNC3",  "1,2,3,4-tetrahydrobenzo[h]isoquinoline"),
    # 1,2,3,4-tetrahydrobenzo[f]isoquinoline: C1-C2-N(3)-C4-C(4a)-C5-C6-C7-C8-C(8a)-C9-C10-C(10a)-C(10b)
    ("CC1CNCc2ccc3ccccc3c21",    "1-methyl-1,2,3,4-tetrahydrobenzo[f]isoquinoline"),
    ("CC1Cc2c(ccc3ccccc23)CN1",  "2-methyl-1,2,3,4-tetrahydrobenzo[f]isoquinoline"),
    ("CN1CCc2c(ccc3ccccc23)C1",  "3-methyl-1,2,3,4-tetrahydrobenzo[f]isoquinoline"),
    ("CC1NCCc2c1ccc1ccccc21",    "4-methyl-1,2,3,4-tetrahydrobenzo[f]isoquinoline"),
    ("Cc1cc2ccccc2c2c1CNCC2",    "5-methyl-1,2,3,4-tetrahydrobenzo[f]isoquinoline"),
    ("Cc1cc2c(c3ccccc13)CCNC2",  "6-methyl-1,2,3,4-tetrahydrobenzo[f]isoquinoline"),
    ("Cc1cccc2c3c(ccc12)CNCC3",  "7-methyl-1,2,3,4-tetrahydrobenzo[f]isoquinoline"),
    ("Cc1ccc2c3c(ccc2c1)CNCC3",  "8-methyl-1,2,3,4-tetrahydrobenzo[f]isoquinoline"),
    ("Cc1ccc2ccc3c(c2c1)CCNC3",  "9-methyl-1,2,3,4-tetrahydrobenzo[f]isoquinoline"),
    ("Cc1cccc2ccc3c(c12)CCNC3",  "10-methyl-1,2,3,4-tetrahydrobenzo[f]isoquinoline"),
    # 1,2,3,4-tetrahydrobenzo[g]isoquinoline: C1-C2-C3-C4-C(4a)-C5-C(5a)-C6-C7-C8-C9-C(9a)-C10-C(10a); N at 2
    ("CC1NCCc2cc3ccccc3cc21",    "1-methyl-1,2,3,4-tetrahydrobenzo[g]isoquinoline"),
    ("CN1CCc2cc3ccccc3cc2C1",    "2-methyl-1,2,3,4-tetrahydrobenzo[g]isoquinoline"),
    ("CC1Cc2cc3ccccc3cc2CN1",    "3-methyl-1,2,3,4-tetrahydrobenzo[g]isoquinoline"),
    ("CC1CNCc2cc3ccccc3cc21",    "4-methyl-1,2,3,4-tetrahydrobenzo[g]isoquinoline"),
    ("Cc1c2c(cc3ccccc13)CNCC2",  "5-methyl-1,2,3,4-tetrahydrobenzo[g]isoquinoline"),
    ("Cc1cccc2cc3c(cc12)CCNC3",  "6-methyl-1,2,3,4-tetrahydrobenzo[g]isoquinoline"),
    ("Cc1ccc2cc3c(cc2c1)CCNC3",  "7-methyl-1,2,3,4-tetrahydrobenzo[g]isoquinoline"),
    ("Cc1ccc2cc3c(cc2c1)CNCC3",  "8-methyl-1,2,3,4-tetrahydrobenzo[g]isoquinoline"),
    ("Cc1cccc2cc3c(cc12)CNCC3",  "9-methyl-1,2,3,4-tetrahydrobenzo[g]isoquinoline"),
    ("Cc1c2c(cc3ccccc13)CCNC2",  "10-methyl-1,2,3,4-tetrahydrobenzo[g]isoquinoline"),
    # 1,2,3,4-tetrahydrobenzo[h]isoquinoline: C1-C2-C3-C4-C(4a)-C5-C6-C7-C8-C(8a)-C9-C10-C(10a)-C(10b); N at 3
    ("CC1NCCc2ccc3ccccc3c21",    "1-methyl-1,2,3,4-tetrahydrobenzo[h]isoquinoline"),
    ("CN1CCc2ccc3ccccc3c2C1",    "2-methyl-1,2,3,4-tetrahydrobenzo[h]isoquinoline"),
    ("CC1Cc2ccc3ccccc3c2CN1",    "3-methyl-1,2,3,4-tetrahydrobenzo[h]isoquinoline"),
    ("CC1CNCc2c1ccc1ccccc21",    "4-methyl-1,2,3,4-tetrahydrobenzo[h]isoquinoline"),
    ("Cc1cc2ccccc2c2c1CCNC2",    "5-methyl-1,2,3,4-tetrahydrobenzo[h]isoquinoline"),
    ("Cc1cc2c(c3ccccc13)CNCC2",  "6-methyl-1,2,3,4-tetrahydrobenzo[h]isoquinoline"),
    ("Cc1cccc2c3c(ccc12)CCNC3",  "7-methyl-1,2,3,4-tetrahydrobenzo[h]isoquinoline"),
    ("Cc1ccc2c3c(ccc2c1)CCNC3",  "8-methyl-1,2,3,4-tetrahydrobenzo[h]isoquinoline"),
    ("Cc1ccc2ccc3c(c2c1)CNCC3",  "9-methyl-1,2,3,4-tetrahydrobenzo[h]isoquinoline"),
    ("Cc1cccc2ccc3c(c12)CNCC3",  "10-methyl-1,2,3,4-tetrahydrobenzo[h]isoquinoline"),
])
def test_phase641_tetrahydrobenzo_isoquinolines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
