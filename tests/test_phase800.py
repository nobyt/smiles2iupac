"""Phase 800: imidazo[1,2-b]pyridazine α-ol/thiol → tautomers.

- C2 → 2(2H)-one/thione; C3 → 3(2H)-one/thione; C6 → 6(4H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imidazo[1,2-b]pyridazine C2-OH/SH
    ("Oc1cn2ncccc2n1",   "imidazo[1,2-b]pyridazin-2(2H)-one"),
    ("Sc1cn2ncccc2n1",   "imidazo[1,2-b]pyridazin-2(2H)-thione"),
    # imidazo[1,2-b]pyridazine C3-OH/SH
    ("Oc1cnc2cccnn12",   "imidazo[1,2-b]pyridazin-3(2H)-one"),
    ("Sc1cnc2cccnn12",   "imidazo[1,2-b]pyridazin-3(2H)-thione"),
    # imidazo[1,2-b]pyridazine C6-OH/SH
    ("Oc1ccc2nccn2n1",   "imidazo[1,2-b]pyridazin-6(4H)-one"),
    ("Sc1ccc2nccn2n1",   "imidazo[1,2-b]pyridazin-6(4H)-thione"),
    # Regression: parent ring unchanged
    ("c1cnn2ccnc2c1",    "imidazo[1,2-b]pyridazine"),
])
def test_phase800_imidazo_pyridazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
