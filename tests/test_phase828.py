"""Phase 828: pyrido[3,4-e]pyrimidine α-ol/thiol → tautomers.

- C2 → 2(1H)-one/thione; C4 → 4(1H)-one/thione
- C5 → 5(1H)-one/thione; C7 → 7(1H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1ncc2cnccc2n1",   "pyrido[3,4-e]pyrimidin-2(1H)-one"),
    ("Sc1ncc2cnccc2n1",   "pyrido[3,4-e]pyrimidin-2(1H)-thione"),
    ("Oc1ncnc2ccncc12",   "pyrido[3,4-e]pyrimidin-4(1H)-one"),
    ("Sc1ncnc2ccncc12",   "pyrido[3,4-e]pyrimidin-4(1H)-thione"),
    ("Oc1nccc2ncncc12",   "pyrido[3,4-e]pyrimidin-5(1H)-one"),
    ("Sc1nccc2ncncc12",   "pyrido[3,4-e]pyrimidin-5(1H)-thione"),
    ("Oc1cc2ncncc2cn1",   "pyrido[3,4-e]pyrimidin-7(1H)-one"),
    ("Sc1cc2ncncc2cn1",   "pyrido[3,4-e]pyrimidin-7(1H)-thione"),
    ("c1cc2ncncc2cn1",    "pyrido[3,4-e]pyrimidine"),
])
def test_phase828(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
