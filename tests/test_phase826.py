"""Phase 826: pyrido[3,4-d]pyrimidine α-ol/thiol → tautomers.

- C2 → 2(1H)-one/thione; C4 → 4(1H)-one/thione
- C6 → 6(2H)-one/thione; C8 → 8(2H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Oc1ncc2ccncc2n1",   "pyrido[3,4-d]pyrimidin-2(1H)-one"),
    ("Sc1ncc2ccncc2n1",   "pyrido[3,4-d]pyrimidin-2(1H)-thione"),
    ("Oc1ncnc2cnccc12",   "pyrido[3,4-d]pyrimidin-4(1H)-one"),
    ("Sc1ncnc2cnccc12",   "pyrido[3,4-d]pyrimidin-4(1H)-thione"),
    ("Oc1cc2cncnc2cn1",   "pyrido[3,4-d]pyrimidin-6(2H)-one"),
    ("Sc1cc2cncnc2cn1",   "pyrido[3,4-d]pyrimidin-6(2H)-thione"),
    ("Oc1nccc2cncnc12",   "pyrido[3,4-d]pyrimidin-8(2H)-one"),
    ("Sc1nccc2cncnc12",   "pyrido[3,4-d]pyrimidin-8(2H)-thione"),
    ("c1cc2cncnc2cn1",    "pyrido[3,4-d]pyrimidine"),
])
def test_phase826(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
