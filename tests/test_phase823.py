"""Phase 823: pyrido[2,3-d]pyridazine α-ol/thiol → tautomers.

- C2 → 2(1H)-one/thione; C5 → 5(1H)-one/thione; C8 → 8(2H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1ccc2cnncc2n1",   "pyrido[2,3-d]pyridazin-2(1H)-one"),
    ("Sc1ccc2cnncc2n1",   "pyrido[2,3-d]pyridazin-2(1H)-thione"),
    ("Oc1nncc2ncccc12",   "pyrido[2,3-d]pyridazin-5(1H)-one"),
    ("Sc1nncc2ncccc12",   "pyrido[2,3-d]pyridazin-5(1H)-thione"),
    ("Oc1nncc2cccnc12",   "pyrido[2,3-d]pyridazin-8(2H)-one"),
    ("Sc1nncc2cccnc12",   "pyrido[2,3-d]pyridazin-8(2H)-thione"),
    ("c1cnc2cnncc2c1",    "pyrido[2,3-d]pyridazine"),
])
def test_phase823(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
