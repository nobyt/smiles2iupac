"""Phase 777: imidazo[1,5-a]-, pyrazolo[1,5-a]-, [1,2,4]triazolo[4,3-a]-pyridine α-ol/thiol → tautomers.

imidazo[1,5-a]pyridine (N2, N4):
- C1 (alpha to N2) → 1(2H)-one/thione
- C3 (alpha to N2 and N4) → 3(2H)-one/thione  [H on lower-locant N2]
- C5 (alpha to N4 junction) → 5(4H)-one/thione

pyrazolo[1,5-a]pyridine (N1, N5 junction):
- C2 (alpha to N5 junction? wait: actually adjacent in 5-membered ring) → 2(5H)-one/thione

[1,2,4]triazolo[4,3-a]pyridine (N1, N2, N4 junction):
- C3 (alpha to N4 junction) → 3(4H)-one/thione
- C5 (alpha to N4 junction, in 6-membered ring) → 5(4H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imidazo[1,5-a]pyridine C1-OH/SH (alpha to N2)
    ("Oc1ncn2ccccc12",      "imidazo[1,5-a]pyridin-1(2H)-one"),
    ("Sc1ncn2ccccc12",      "imidazo[1,5-a]pyridin-1(2H)-thione"),
    # imidazo[1,5-a]pyridine C3-OH/SH (alpha to N2 and N4; H on N2)
    ("Oc1ncc2ccccn12",      "imidazo[1,5-a]pyridin-3(2H)-one"),
    ("Sc1ncc2ccccn12",      "imidazo[1,5-a]pyridin-3(2H)-thione"),
    # imidazo[1,5-a]pyridine C5-OH/SH (alpha to N4 junction)
    ("Oc1cccc2cncn12",      "imidazo[1,5-a]pyridin-5(4H)-one"),
    ("Sc1cccc2cncn12",      "imidazo[1,5-a]pyridin-5(4H)-thione"),
    # pyrazolo[1,5-a]pyridine C2-OH/SH
    ("Oc1cc2ccccn2n1",      "pyrazolo[1,5-a]pyridin-2(5H)-one"),
    ("Sc1cc2ccccn2n1",      "pyrazolo[1,5-a]pyridin-2(5H)-thione"),
    # [1,2,4]triazolo[4,3-a]pyridine C3-OH/SH (alpha to N4 junction)
    ("Oc1nnc2ccccn12",      "[1,2,4]triazolo[4,3-a]pyridin-3(4H)-one"),
    ("Sc1nnc2ccccn12",      "[1,2,4]triazolo[4,3-a]pyridin-3(4H)-thione"),
    # [1,2,4]triazolo[4,3-a]pyridine C5-OH/SH (alpha to N4 junction, 6-membered ring)
    ("Oc1cccc2nncn12",      "[1,2,4]triazolo[4,3-a]pyridin-5(4H)-one"),
    ("Sc1cccc2nncn12",      "[1,2,4]triazolo[4,3-a]pyridin-5(4H)-thione"),
    # Regressions: parent rings unchanged
    ("c1ccn2cncc2c1",       "imidazo[1,5-a]pyridine"),
    ("c1ccn2nccc2c1",       "pyrazolo[1,5-a]pyridine"),
    ("c1ccn2cnnc2c1",       "[1,2,4]triazolo[4,3-a]pyridine"),
])
def test_phase777_fused_5_6_az_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
