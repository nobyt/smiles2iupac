"""Phase 813: [1,2,4]triazolo[4,3-b][1,2,4]triazine α-ol/thiol → tautomers.

- C3 → 3(2H)-one/thione; C6 → 6(4H)-one/thione; C7 → 7(1H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1nnc2nccnn12",   "[1,2,4]triazolo[4,3-b][1,2,4]triazin-3(2H)-one"),
    ("Sc1nnc2nccnn12",   "[1,2,4]triazolo[4,3-b][1,2,4]triazin-3(2H)-thione"),
    ("Oc1cnc2nncn2n1",   "[1,2,4]triazolo[4,3-b][1,2,4]triazin-6(4H)-one"),
    ("Sc1cnc2nncn2n1",   "[1,2,4]triazolo[4,3-b][1,2,4]triazin-6(4H)-thione"),
    ("Oc1cnn2cnnc2n1",   "[1,2,4]triazolo[4,3-b][1,2,4]triazin-7(1H)-one"),
    ("Sc1cnn2cnnc2n1",   "[1,2,4]triazolo[4,3-b][1,2,4]triazin-7(1H)-thione"),
    ("c1cnn2cnnc2n1",    "[1,2,4]triazolo[4,3-b][1,2,4]triazine"),
])
def test_phase813(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
