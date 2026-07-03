"""Phase 789: 1H-pyrazolo[4,5-e]pyrazine α-ol/thiol → tautomers.

- C3 (alpha to N2) → 3(2H)-one/thione
- C5 (alpha to junction N4) → 5(4H)-one/thione
- C6 (alpha to junction N4) → 6(5H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # C3-OH/SH (alpha to N2)
    ("Oc1n[nH]c2nccnc12",   "1H-pyrazolo[4,5-e]pyrazin-3(2H)-one"),
    ("Sc1n[nH]c2nccnc12",   "1H-pyrazolo[4,5-e]pyrazin-3(2H)-thione"),
    # C5-OH/SH (alpha to junction N4)
    ("Oc1cnc2[nH]ncc2n1",   "1H-pyrazolo[4,5-e]pyrazin-5(4H)-one"),
    ("Sc1cnc2[nH]ncc2n1",   "1H-pyrazolo[4,5-e]pyrazin-5(4H)-thione"),
    # C6-OH/SH (alpha to junction N5)
    ("Oc1cnc2cn[nH]c2n1",   "1H-pyrazolo[4,5-e]pyrazin-6(5H)-one"),
    ("Sc1cnc2cn[nH]c2n1",   "1H-pyrazolo[4,5-e]pyrazin-6(5H)-thione"),
    # Regression: parent unchanged
    ("c1n[nH]c2nccnc12",    "1H-pyrazolo[4,5-e]pyrazine"),
])
def test_phase789_pyrazolo_pyrazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
