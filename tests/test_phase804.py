"""Phase 804: pyrazolo[1,5-a]pyrazine α-ol/thiol → tautomers.

- C2 → 2(1H)-one/thione; C4 → 4(5H)-one/thione
- C6 → 6(6H)-one/thione; C7 → 7(6H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrazolo[1,5-a]pyrazine C2-OH/SH
    ("Oc1cc2cnccn2n1",   "pyrazolo[1,5-a]pyrazin-2(1H)-one"),
    ("Sc1cc2cnccn2n1",   "pyrazolo[1,5-a]pyrazin-2(1H)-thione"),
    # pyrazolo[1,5-a]pyrazine C4-OH/SH
    ("Oc1nccn2nccc12",   "pyrazolo[1,5-a]pyrazin-4(5H)-one"),
    ("Sc1nccn2nccc12",   "pyrazolo[1,5-a]pyrazin-4(5H)-thione"),
    # pyrazolo[1,5-a]pyrazine C6-OH/SH
    ("Oc1cn2nccc2cn1",   "pyrazolo[1,5-a]pyrazin-6(6H)-one"),
    ("Sc1cn2nccc2cn1",   "pyrazolo[1,5-a]pyrazin-6(6H)-thione"),
    # pyrazolo[1,5-a]pyrazine C7-OH/SH
    ("Oc1cncc2ccnn12",   "pyrazolo[1,5-a]pyrazin-7(6H)-one"),
    ("Sc1cncc2ccnn12",   "pyrazolo[1,5-a]pyrazin-7(6H)-thione"),
    # Regression: parent ring unchanged
    ("c1cn2nccc2cn1",    "pyrazolo[1,5-a]pyrazine"),
])
def test_phase804_pyrazolo_pyrazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
