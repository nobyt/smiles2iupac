"""Phase 792: 1H-pyrrolo[2,3-c]pyridine, 1H-pyrrolo[3,2-c]pyridine α-ol/thiol → tautomers.

- 1H-pyrrolo[2,3-c]pyridine C2 → 2(3H)-one/thione; C5 → 5(4H)-one/thione; C7 → 7(6H)-one/thione
- 1H-pyrrolo[3,2-c]pyridine C2 → 2(3H)-one/thione; C4 → 4(3H)-one/thione; C6 → 6(5H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-pyrrolo[2,3-c]pyridine C2-OH/SH
    ("Oc1cc2ccncc2[nH]1",  "1H-pyrrolo[2,3-c]pyridin-2(3H)-one"),
    ("Sc1cc2ccncc2[nH]1",  "1H-pyrrolo[2,3-c]pyridin-2(3H)-thione"),
    # 1H-pyrrolo[2,3-c]pyridine C5-OH/SH
    ("Oc1cc2cc[nH]c2cn1",  "1H-pyrrolo[2,3-c]pyridin-5(4H)-one"),
    ("Sc1cc2cc[nH]c2cn1",  "1H-pyrrolo[2,3-c]pyridin-5(4H)-thione"),
    # 1H-pyrrolo[2,3-c]pyridine C7-OH/SH
    ("Oc1nccc2cc[nH]c12",  "1H-pyrrolo[2,3-c]pyridin-7(6H)-one"),
    ("Sc1nccc2cc[nH]c12",  "1H-pyrrolo[2,3-c]pyridin-7(6H)-thione"),
    # 1H-pyrrolo[3,2-c]pyridine C2-OH/SH
    ("Oc1cc2cnccc2[nH]1",  "1H-pyrrolo[3,2-c]pyridin-2(3H)-one"),
    ("Sc1cc2cnccc2[nH]1",  "1H-pyrrolo[3,2-c]pyridin-2(3H)-thione"),
    # 1H-pyrrolo[3,2-c]pyridine C4-OH/SH
    ("Oc1nccc2[nH]ccc12",  "1H-pyrrolo[3,2-c]pyridin-4(3H)-one"),
    ("Sc1nccc2[nH]ccc12",  "1H-pyrrolo[3,2-c]pyridin-4(3H)-thione"),
    # 1H-pyrrolo[3,2-c]pyridine C6-OH/SH
    ("Oc1cc2[nH]ccc2cn1",  "1H-pyrrolo[3,2-c]pyridin-6(5H)-one"),
    ("Sc1cc2[nH]ccc2cn1",  "1H-pyrrolo[3,2-c]pyridin-6(5H)-thione"),
    # Regressions: parent rings unchanged
    ("c1cncc2[nH]ccc12",   "1H-pyrrolo[2,3-c]pyridine"),
    ("c1cncc2cc[nH]c12",   "1H-pyrrolo[3,2-c]pyridine"),
])
def test_phase792_pyrrolo_pyridine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
