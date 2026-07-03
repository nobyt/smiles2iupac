"""Phase 816: imidazo[1,5-a]pyrazine α-ol/thiol → tautomers.

- C1 → 1(5H)-one/thione; C3 → 3(2H)-one/thione; C5 → 5(1H)-one/thione
- C6 → 6(4H)-one/thione; C8 → 8(5H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1ncn2ccncc12",   "imidazo[1,5-a]pyrazin-1(5H)-one"),
    ("Sc1ncn2ccncc12",   "imidazo[1,5-a]pyrazin-1(5H)-thione"),
    ("Oc1ncc2cnccn12",   "imidazo[1,5-a]pyrazin-3(2H)-one"),
    ("Sc1ncc2cnccn12",   "imidazo[1,5-a]pyrazin-3(2H)-thione"),
    ("Oc1cncc2cncn12",   "imidazo[1,5-a]pyrazin-5(1H)-one"),
    ("Sc1cncc2cncn12",   "imidazo[1,5-a]pyrazin-5(1H)-thione"),
    ("Oc1cn2cncc2cn1",   "imidazo[1,5-a]pyrazin-6(4H)-one"),
    ("Sc1cn2cncc2cn1",   "imidazo[1,5-a]pyrazin-6(4H)-thione"),
    ("Oc1nccn2cncc12",   "imidazo[1,5-a]pyrazin-8(5H)-one"),
    ("Sc1nccn2cncc12",   "imidazo[1,5-a]pyrazin-8(5H)-thione"),
    ("c1cn2cncc2cn1",    "imidazo[1,5-a]pyrazine"),
])
def test_phase816(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
