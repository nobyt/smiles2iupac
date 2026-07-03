"""Phase 827: pyrido[3,4-e]pyridazine α-ol/thiol → tautomers.

- C3 → 3(2H)-one/thione; C5 → 5(1H)-one/thione; C7 → 7(1H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1cc2cnccc2nn1",   "pyrido[3,4-e]pyridazin-3(2H)-one"),
    ("Sc1cc2cnccc2nn1",   "pyrido[3,4-e]pyridazin-3(2H)-thione"),
    ("Oc1nccc2nnccc12",   "pyrido[3,4-e]pyridazin-5(1H)-one"),
    ("Sc1nccc2nnccc12",   "pyrido[3,4-e]pyridazin-5(1H)-thione"),
    ("Oc1cc2nnccc2cn1",   "pyrido[3,4-e]pyridazin-7(1H)-one"),
    ("Sc1cc2nnccc2cn1",   "pyrido[3,4-e]pyridazin-7(1H)-thione"),
    ("c1cc2nnccc2cn1",    "pyrido[3,4-e]pyridazine"),
])
def test_phase827(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
