"""Phase 810: [1,2,3]triazolo[1,5-a]pyrimidine α-ol/thiol → tautomers.

- C3 → 3(5H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1nnn2cccnc12",   "[1,2,3]triazolo[1,5-a]pyrimidin-3(5H)-one"),
    ("Sc1nnn2cccnc12",   "[1,2,3]triazolo[1,5-a]pyrimidin-3(5H)-thione"),
    ("c1cnc2cnnn2c1",    "[1,2,3]triazolo[1,5-a]pyrimidine"),
])
def test_phase810(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
