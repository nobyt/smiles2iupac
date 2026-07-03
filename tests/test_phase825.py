"""Phase 825: pyrido[3,4-c]pyridazine α-ol/thiol → tautomers.

- C3 → 3(2H)-one/thione; C6 → 6(2H)-one/thione; C8 → 8(2H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1cc2ccncc2nn1",   "pyrido[3,4-c]pyridazin-3(2H)-one"),
    ("Sc1cc2ccncc2nn1",   "pyrido[3,4-c]pyridazin-3(2H)-thione"),
    ("Oc1cc2ccnnc2cn1",   "pyrido[3,4-c]pyridazin-6(2H)-one"),
    ("Sc1cc2ccnnc2cn1",   "pyrido[3,4-c]pyridazin-6(2H)-thione"),
    ("Oc1nccc2ccnnc12",   "pyrido[3,4-c]pyridazin-8(2H)-one"),
    ("Sc1nccc2ccnnc12",   "pyrido[3,4-c]pyridazin-8(2H)-thione"),
    ("c1cc2ccnnc2cn1",    "pyrido[3,4-c]pyridazine"),
])
def test_phase825(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
