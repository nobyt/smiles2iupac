"""Phase 515: 1,2-oxazinane and 1,2-thiazinane naming
(IUPAC 2013 P-22.2: Hantzsch-Widman names for saturated 6-membered rings with
adjacent O+N or S+N heteroatoms).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2-oxazinane (O at 1, N at 2, 6-membered saturated)
    ("O1NCCCC1",  "1,2-oxazinane"),
    # 1,2-thiazinane (S at 1, N at 2, 6-membered saturated)
    ("S1NCCCC1",  "1,2-thiazinane"),
    # Regression: other oxazinane/thiazinane regioisomers unaffected
    ("O1CCNCC1",  "morpholine"),         # 1,4-oxazinane (retained name)
    ("S1CCNCC1",  "1,4-thiazinane"),
    ("S1CCCNC1",  "1,3-thiazinane"),
])
def test_phase515_1_2_oxazinane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
