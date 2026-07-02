"""Phase 778: imidazo[1,2-a]pyridine α-ol/thiol → tautomers (IUPAC 2013).

imidazo[1,2-a]pyridine has two N atoms: non-junction N and junction N (N4).
Alpha C positions:
- C2 (in 5-membered ring, alpha to non-junction N) → 2(3H)-one/thione
- C3 (in 5-membered ring, alpha to junction N and non-junction N; H on lower locant) → 3(2H)-one/thione
- C5 (in 6-membered ring, alpha to junction N4) → 5(4H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # C2-OH/SH (alpha to non-junction N, H on adjacent N3)
    ("Oc1cn2ccccc2n1",      "imidazo[1,2-a]pyridin-2(3H)-one"),
    ("Sc1cn2ccccc2n1",      "imidazo[1,2-a]pyridin-2(3H)-thione"),
    # C3-OH/SH (alpha to junction N4; H on N2 = lower locant)
    ("Oc1cnc2ccccn12",      "imidazo[1,2-a]pyridin-3(2H)-one"),
    ("Sc1cnc2ccccn12",      "imidazo[1,2-a]pyridin-3(2H)-thione"),
    # C5-OH/SH (alpha to junction N4 in 6-membered ring)
    ("Oc1cccc2nccn12",      "imidazo[1,2-a]pyridin-5(4H)-one"),
    ("Sc1cccc2nccn12",      "imidazo[1,2-a]pyridin-5(4H)-thione"),
    # Regression: parent ring unchanged
    ("c1ccn2ccnc2c1",       "imidazo[1,2-a]pyridine"),
])
def test_phase778_imidazo_12a_pyridine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
