"""Phase 803: [1,2,4]triazolo[4,3-a]pyrimidine α-ol/thiol → tautomers.

- C3 → 3(2H)-one/thione; C5 → 5(6H)-one/thione; C7 → 7(8H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # [1,2,4]triazolo[4,3-a]pyrimidine C3-OH/SH
    ("Oc1nnc2ncccn12",   "[1,2,4]triazolo[4,3-a]pyrimidin-3(2H)-one"),
    ("Sc1nnc2ncccn12",   "[1,2,4]triazolo[4,3-a]pyrimidin-3(2H)-thione"),
    # [1,2,4]triazolo[4,3-a]pyrimidine C5-OH/SH
    ("Oc1ccnc2nncn12",   "[1,2,4]triazolo[4,3-a]pyrimidin-5(6H)-one"),
    ("Sc1ccnc2nncn12",   "[1,2,4]triazolo[4,3-a]pyrimidin-5(6H)-thione"),
    # [1,2,4]triazolo[4,3-a]pyrimidine C7-OH/SH
    ("Oc1ccn2cnnc2n1",   "[1,2,4]triazolo[4,3-a]pyrimidin-7(8H)-one"),
    ("Sc1ccn2cnnc2n1",   "[1,2,4]triazolo[4,3-a]pyrimidin-7(8H)-thione"),
    # Regression: parent ring unchanged
    ("c1cnc2nncn2c1",    "[1,2,4]triazolo[4,3-a]pyrimidine"),
])
def test_phase803_triazolo_pyrimidine_43a_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
