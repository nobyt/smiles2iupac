"""Phase 797: imidazo[1,2-c]pyrimidine α-ol/thiol → tautomers.

- C2 → 2(1H)-one/thione (adjacent to N1)
- C3 → 3(4H)-one/thione (adjacent to N4)
- C5 → 5(4H)-one/thione (adjacent to N4, lower locant than N6)
- C7 → 7(6H)-one/thione (adjacent to N6)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imidazo[1,2-c]pyrimidine C2-OH/SH
    ("Oc1cn2cnccc2n1",   "imidazo[1,2-c]pyrimidin-2(1H)-one"),
    ("Sc1cn2cnccc2n1",   "imidazo[1,2-c]pyrimidin-2(1H)-thione"),
    # imidazo[1,2-c]pyrimidine C3-OH/SH
    ("Oc1cnc2ccncn12",   "imidazo[1,2-c]pyrimidin-3(4H)-one"),
    ("Sc1cnc2ccncn12",   "imidazo[1,2-c]pyrimidin-3(4H)-thione"),
    # imidazo[1,2-c]pyrimidine C5-OH/SH
    ("Oc1nccc2nccn12",   "imidazo[1,2-c]pyrimidin-5(4H)-one"),
    ("Sc1nccc2nccn12",   "imidazo[1,2-c]pyrimidin-5(4H)-thione"),
    # imidazo[1,2-c]pyrimidine C7-OH/SH
    ("Oc1cc2nccn2cn1",   "imidazo[1,2-c]pyrimidin-7(6H)-one"),
    ("Sc1cc2nccn2cn1",   "imidazo[1,2-c]pyrimidin-7(6H)-thione"),
    # Regression: parent ring unchanged
    ("c1cnc2ccncn12",    "imidazo[1,2-c]pyrimidine"),
])
def test_phase797_imidazo_pyrimidine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
