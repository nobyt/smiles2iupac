"""Phase 808: [1,2,4]triazolo[4,3-a]pyrazine α-ol/thiol → tautomers.

- C3 → 3(2H)-one/thione; C5 → 5(1H)-one/thione
- C6 → 6(4H)-one/thione; C8 → 8(5H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # [1,2,4]triazolo[4,3-a]pyrazine C3-OH/SH
    ("Oc1nnc2cnccn12",   "[1,2,4]triazolo[4,3-a]pyrazin-3(2H)-one"),
    ("Sc1nnc2cnccn12",   "[1,2,4]triazolo[4,3-a]pyrazin-3(2H)-thione"),
    # [1,2,4]triazolo[4,3-a]pyrazine C5-OH/SH
    ("Oc1cncc2nncn12",   "[1,2,4]triazolo[4,3-a]pyrazin-5(1H)-one"),
    ("Sc1cncc2nncn12",   "[1,2,4]triazolo[4,3-a]pyrazin-5(1H)-thione"),
    # [1,2,4]triazolo[4,3-a]pyrazine C6-OH/SH
    ("Oc1cn2cnnc2cn1",   "[1,2,4]triazolo[4,3-a]pyrazin-6(4H)-one"),
    ("Sc1cn2cnnc2cn1",   "[1,2,4]triazolo[4,3-a]pyrazin-6(4H)-thione"),
    # [1,2,4]triazolo[4,3-a]pyrazine C8-OH/SH
    ("Oc1nccn2cnnc12",   "[1,2,4]triazolo[4,3-a]pyrazin-8(5H)-one"),
    ("Sc1nccn2cnnc12",   "[1,2,4]triazolo[4,3-a]pyrazin-8(5H)-thione"),
    # Regression: parent ring unchanged
    ("c1cn2cnnc2cn1",    "[1,2,4]triazolo[4,3-a]pyrazine"),
])
def test_phase808_triazolo_pyrazine_43a_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
