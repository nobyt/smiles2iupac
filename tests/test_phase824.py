"""Phase 824: pyrido[2,3-e]pyridazine α-ol/thiol → tautomers.

- C3 → 3(2H)-one/thione; C6 → 6(2H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1cc2ncccc2nn1",   "pyrido[2,3-e]pyridazin-3(2H)-one"),
    ("Sc1cc2ncccc2nn1",   "pyrido[2,3-e]pyridazin-3(2H)-thione"),
    ("Oc1ccc2nnccc2n1",   "pyrido[2,3-e]pyridazin-6(2H)-one"),
    ("Sc1ccc2nnccc2n1",   "pyrido[2,3-e]pyridazin-6(2H)-thione"),
    ("c1cnc2ccnnc2c1",    "pyrido[2,3-e]pyridazine"),
])
def test_phase824(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
