"""Phase 779: imidazo[1,2-a]pyrimidine α-ol/thiol → tautomers (IUPAC 2013).

imidazo[1,2-a]pyrimidine alpha C positions:
- C2 (in 5-membered ring, alpha to N3) → 2(3H)-one/thione
- C5 (in 6-membered ring, alpha to junction N4) → 5(4H)-one/thione
- C7 (in 6-membered ring, alpha to junction N4) → 7(4H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # C2-OH/SH (alpha to N3)
    ("Oc1cn2cccnc2n1",      "imidazo[1,2-a]pyrimidin-2(3H)-one"),
    ("Sc1cn2cccnc2n1",      "imidazo[1,2-a]pyrimidin-2(3H)-thione"),
    # C5-OH/SH (alpha to junction N4)
    ("Oc1ccnc2nccn12",      "imidazo[1,2-a]pyrimidin-5(4H)-one"),
    ("Sc1ccnc2nccn12",      "imidazo[1,2-a]pyrimidin-5(4H)-thione"),
    # C7-OH/SH (alpha to junction N4)
    ("Oc1ccn2ccnc2n1",      "imidazo[1,2-a]pyrimidin-7(4H)-one"),
    ("Sc1ccn2ccnc2n1",      "imidazo[1,2-a]pyrimidin-7(4H)-thione"),
    # Regression: parent ring unchanged
    ("c1cnc2nccn2c1",       "imidazo[1,2-a]pyrimidine"),
])
def test_phase779_imidazo_12a_pyrimidine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
