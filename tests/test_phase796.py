"""Phase 796: pyrido[2,3-d]pyrimidine, pyrido[3,4-e]pyrazine α-ol/thiol → tautomers.

- pyrido[2,3-d]pyrimidine C2 → 2(1H)-one/thione; C4 → 4(3H)-one/thione; C7 → 7(8H)-one/thione
- pyrido[3,4-e]pyrazine C2 → 2(1H)-one/thione; C3 → 3(4H)-one/thione; C5 → 5(6H)-one/thione; C7 → 7(8H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrido[2,3-d]pyrimidine C2-OH/SH
    ("Oc1ncc2cccnc2n1",   "pyrido[2,3-d]pyrimidin-2(1H)-one"),
    ("Sc1ncc2cccnc2n1",   "pyrido[2,3-d]pyrimidin-2(1H)-thione"),
    # pyrido[2,3-d]pyrimidine C4-OH/SH
    ("Oc1ncnc2ncccc12",   "pyrido[2,3-d]pyrimidin-4(3H)-one"),
    ("Sc1ncnc2ncccc12",   "pyrido[2,3-d]pyrimidin-4(3H)-thione"),
    # pyrido[2,3-d]pyrimidine C7-OH/SH
    ("Oc1ccc2cncnc2n1",   "pyrido[2,3-d]pyrimidin-7(8H)-one"),
    ("Sc1ccc2cncnc2n1",   "pyrido[2,3-d]pyrimidin-7(8H)-thione"),
    # pyrido[3,4-e]pyrazine C2-OH/SH
    ("Oc1cnc2cnccc2n1",   "pyrido[3,4-e]pyrazin-2(1H)-one"),
    ("Sc1cnc2cnccc2n1",   "pyrido[3,4-e]pyrazin-2(1H)-thione"),
    # pyrido[3,4-e]pyrazine C3-OH/SH
    ("Oc1cnc2ccncc2n1",   "pyrido[3,4-e]pyrazin-3(4H)-one"),
    ("Sc1cnc2ccncc2n1",   "pyrido[3,4-e]pyrazin-3(4H)-thione"),
    # pyrido[3,4-e]pyrazine C5-OH/SH
    ("Oc1nccc2nccnc12",   "pyrido[3,4-e]pyrazin-5(6H)-one"),
    ("Sc1nccc2nccnc12",   "pyrido[3,4-e]pyrazin-5(6H)-thione"),
    # pyrido[3,4-e]pyrazine C7-OH/SH
    ("Oc1cc2nccnc2cn1",   "pyrido[3,4-e]pyrazin-7(8H)-one"),
    ("Sc1cc2nccnc2cn1",   "pyrido[3,4-e]pyrazin-7(8H)-thione"),
    # Regressions: parent rings unchanged
    ("c1ccnc2ncncc12",    "pyrido[2,3-d]pyrimidine"),
    ("c1cncc2nccnc12",    "pyrido[3,4-e]pyrazine"),
])
def test_phase796_pyrido_pyrimidine_pyrazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
