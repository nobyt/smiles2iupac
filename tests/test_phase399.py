"""Phase 399: Thiolactam indicated-hydrogen naming (IUPAC 2013 P-31.1.7).

Two fixes:
1. heterocycle_handler.py thiolactam block: apply same multi-N IUPAC rotation
   algorithm as lactam block, giving {ring}-{loc}({nh}H)-thione format.
2. group_namers._name_thiourea_if_match: skip in-ring C atoms so thiolactam
   C=S (in ring) is not misidentified as thiourea.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyridin-2(1H)-thione
    ("S=C1NC=CC=C1",    "pyridin-2(1H)-thione"),
    # pyrimidin-2(1H)-thione
    ("S=C1NC=CC=N1",    "pyrimidin-2(1H)-thione"),
    # pyrimidin-4(3H)-thione
    ("S=C1NC=NC=C1",    "pyrimidin-4(3H)-thione"),
    # pyrazin-2(1H)-thione
    ("S=C1NC=CN=C1",    "pyrazin-2(1H)-thione"),
    # regression: thiourea still works
    ("NC(=S)N",          "thiourea"),
    # regression: saturated lactam unchanged
    ("O=C1CCCCN1",       "piperidin-2-one"),
    # regression: plain parent rings unchanged
    ("c1ccccn1",          "pyridine"),
    ("c1ccncn1",          "pyrimidine"),
])
def test_phase399_thiolactam(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
