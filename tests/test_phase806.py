"""Phase 806: pyrrolo[1,2-a]pyrazine α-ol/thiol → tautomers.

- C1 → 1(1H)-one/thione; C3 → 3(3H)-one/thione
- C4 → 4(1H)-one/thione; C6 → 6(2H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrrolo[1,2-a]pyrazine C1-OH/SH
    ("Oc1nccn2cccc12",   "pyrrolo[1,2-a]pyrazin-1(1H)-one"),
    ("Sc1nccn2cccc12",   "pyrrolo[1,2-a]pyrazin-1(1H)-thione"),
    # pyrrolo[1,2-a]pyrazine C3-OH/SH
    ("Oc1cn2cccc2cn1",   "pyrrolo[1,2-a]pyrazin-3(3H)-one"),
    ("Sc1cn2cccc2cn1",   "pyrrolo[1,2-a]pyrazin-3(3H)-thione"),
    # pyrrolo[1,2-a]pyrazine C4-OH/SH
    ("Oc1cncc2cccn12",   "pyrrolo[1,2-a]pyrazin-4(1H)-one"),
    ("Sc1cncc2cccn12",   "pyrrolo[1,2-a]pyrazin-4(1H)-thione"),
    # pyrrolo[1,2-a]pyrazine C6-OH/SH
    ("Oc1ccc2cnccn12",   "pyrrolo[1,2-a]pyrazin-6(2H)-one"),
    ("Sc1ccc2cnccn12",   "pyrrolo[1,2-a]pyrazin-6(2H)-thione"),
    # Regression: parent ring unchanged
    ("c1cc2cnccn2c1",    "pyrrolo[1,2-a]pyrazine"),
])
def test_phase806_pyrrolo_pyrazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
