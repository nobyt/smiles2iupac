"""Phase 817: imidazo[1,5-b][1,2,4]triazine α-ol/thiol → tautomers.

- C2 → 2(1H)-one/thione; C3 → 3(3H)-one/thione
- C6 → 6(2H)-one/thione; C8 → 8(2H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1cnn2cncc2n1",   "imidazo[1,5-b][1,2,4]triazin-2(1H)-one"),
    ("Sc1cnn2cncc2n1",   "imidazo[1,5-b][1,2,4]triazin-2(1H)-thione"),
    ("Oc1cnc2cncn2n1",   "imidazo[1,5-b][1,2,4]triazin-3(3H)-one"),
    ("Sc1cnc2cncn2n1",   "imidazo[1,5-b][1,2,4]triazin-3(3H)-thione"),
    ("Oc1ncc2nccnn12",   "imidazo[1,5-b][1,2,4]triazin-6(2H)-one"),
    ("Sc1ncc2nccnn12",   "imidazo[1,5-b][1,2,4]triazin-6(2H)-thione"),
    ("Oc1ncn2nccnc12",   "imidazo[1,5-b][1,2,4]triazin-8(2H)-one"),
    ("Sc1ncn2nccnc12",   "imidazo[1,5-b][1,2,4]triazin-8(2H)-thione"),
    ("c1cnn2cncc2n1",    "imidazo[1,5-b][1,2,4]triazine"),
])
def test_phase817(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
