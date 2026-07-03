"""Phase 801: imidazo[1,2-a]pyrazine α-ol/thiol → tautomers.

- C2 → 2(2H)-one/thione; C3 → 3(2H)-one/thione; C5 → 5(6H)-one/thione
- C6 → 6(4H)-one/thione; C8 → 8(7H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imidazo[1,2-a]pyrazine C2-OH/SH
    ("Oc1cn2ccncc2n1",   "imidazo[1,2-a]pyrazin-2(2H)-one"),
    ("Sc1cn2ccncc2n1",   "imidazo[1,2-a]pyrazin-2(2H)-thione"),
    # imidazo[1,2-a]pyrazine C3-OH/SH
    ("Oc1cnc2cnccn12",   "imidazo[1,2-a]pyrazin-3(2H)-one"),
    ("Sc1cnc2cnccn12",   "imidazo[1,2-a]pyrazin-3(2H)-thione"),
    # imidazo[1,2-a]pyrazine C5-OH/SH
    ("Oc1cncc2nccn12",   "imidazo[1,2-a]pyrazin-5(6H)-one"),
    ("Sc1cncc2nccn12",   "imidazo[1,2-a]pyrazin-5(6H)-thione"),
    # imidazo[1,2-a]pyrazine C6-OH/SH
    ("Oc1cn2ccnc2cn1",   "imidazo[1,2-a]pyrazin-6(4H)-one"),
    ("Sc1cn2ccnc2cn1",   "imidazo[1,2-a]pyrazin-6(4H)-thione"),
    # imidazo[1,2-a]pyrazine C8-OH/SH
    ("Oc1nccn2ccnc12",   "imidazo[1,2-a]pyrazin-8(7H)-one"),
    ("Sc1nccn2ccnc12",   "imidazo[1,2-a]pyrazin-8(7H)-thione"),
    # Regression: parent ring unchanged
    ("c1cn2ccnc2cn1",    "imidazo[1,2-a]pyrazine"),
])
def test_phase801_imidazo_pyrazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
