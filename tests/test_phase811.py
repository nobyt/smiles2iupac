"""Phase 811: [1,2,3]triazolo[1,5-b]pyridazine α-ol/thiol → tautomers.

- C3 → 3(5H)-one/thione; C6 → 6(6H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1nnn2ncccc12",   "[1,2,3]triazolo[1,5-b]pyridazin-3(5H)-one"),
    ("Sc1nnn2ncccc12",   "[1,2,3]triazolo[1,5-b]pyridazin-3(5H)-thione"),
    ("Oc1ccc2cnnn2n1",   "[1,2,3]triazolo[1,5-b]pyridazin-6(6H)-one"),
    ("Sc1ccc2cnnn2n1",   "[1,2,3]triazolo[1,5-b]pyridazin-6(6H)-thione"),
    ("c1cnn2nncc2c1",    "[1,2,3]triazolo[1,5-b]pyridazine"),
])
def test_phase811(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
