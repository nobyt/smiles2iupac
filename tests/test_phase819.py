"""Phase 819: imidazo[3,2-b][1,2,4]triazine α-ol/thiol → tautomers.

- C2 → 2(1H)-one/thione; C3 → 3(4H)-one/thione
- C6 → 6(6H)-one/thione; C7 → 7(1H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1cnc2nccn2n1",   "imidazo[3,2-b][1,2,4]triazin-2(1H)-one"),
    ("Sc1cnc2nccn2n1",   "imidazo[3,2-b][1,2,4]triazin-2(1H)-thione"),
    ("Oc1cnn2ccnc2n1",   "imidazo[3,2-b][1,2,4]triazin-3(4H)-one"),
    ("Sc1cnn2ccnc2n1",   "imidazo[3,2-b][1,2,4]triazin-3(4H)-thione"),
    ("Oc1cn2nccnc2n1",   "imidazo[3,2-b][1,2,4]triazin-6(6H)-one"),
    ("Sc1cn2nccnc2n1",   "imidazo[3,2-b][1,2,4]triazin-6(6H)-thione"),
    ("Oc1cnc2nccnn12",   "imidazo[3,2-b][1,2,4]triazin-7(1H)-one"),
    ("Sc1cnc2nccnn12",   "imidazo[3,2-b][1,2,4]triazin-7(1H)-thione"),
    ("c1cnn2ccnc2n1",    "imidazo[3,2-b][1,2,4]triazine"),
])
def test_phase819(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
