"""Phase 809: [1,2,3]triazolo[1,5-a]pyridine α-ol/thiol → tautomers.

- C3 → 3(5H)-one/thione; C7 → 7(1H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # [1,2,3]triazolo[1,5-a]pyridine C3-OH/SH
    ("Oc1nnn2ccccc12",   "[1,2,3]triazolo[1,5-a]pyridin-3(5H)-one"),
    ("Sc1nnn2ccccc12",   "[1,2,3]triazolo[1,5-a]pyridin-3(5H)-thione"),
    # [1,2,3]triazolo[1,5-a]pyridine C7-OH/SH
    ("Oc1cccc2cnnn12",   "[1,2,3]triazolo[1,5-a]pyridin-7(1H)-one"),
    ("Sc1cccc2cnnn12",   "[1,2,3]triazolo[1,5-a]pyridin-7(1H)-thione"),
    # Regression: parent ring unchanged
    ("c1ccn2nncc2c1",    "[1,2,3]triazolo[1,5-a]pyridine"),
])
def test_phase809_triazolo_pyridine_123_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
