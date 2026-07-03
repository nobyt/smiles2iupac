"""Phase 793: 1H-pyrazolo[5,4-c]pyridine, 1H-pyrazolo[4,5-c]pyridine, 1H-imidazo[4,5-c]pyridine α-ol/thiol → tautomers.

- 1H-pyrazolo[5,4-c]pyridine C3 → 3(2H)-one/thione; C5 → 5(4H)-one/thione; C7 → 7(6H)-one/thione
- 1H-pyrazolo[4,5-c]pyridine C3 → 3(2H)-one/thione; C4 → 4(3H)-one/thione; C6 → 6(5H)-one/thione
- 1H-imidazo[4,5-c]pyridine C2 → 2(3H)-one/thione; C4 → 4(3H)-one/thione; C6 → 6(5H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-pyrazolo[5,4-c]pyridine C3-OH/SH
    ("Oc1n[nH]c2cnccc12",  "1H-pyrazolo[5,4-c]pyridin-3(2H)-one"),
    ("Sc1n[nH]c2cnccc12",  "1H-pyrazolo[5,4-c]pyridin-3(2H)-thione"),
    # 1H-pyrazolo[5,4-c]pyridine C5-OH/SH
    ("Oc1cc2cn[nH]c2cn1",  "1H-pyrazolo[5,4-c]pyridin-5(4H)-one"),
    ("Sc1cc2cn[nH]c2cn1",  "1H-pyrazolo[5,4-c]pyridin-5(4H)-thione"),
    # 1H-pyrazolo[5,4-c]pyridine C7-OH/SH
    ("Oc1nccc2cn[nH]c12",  "1H-pyrazolo[5,4-c]pyridin-7(6H)-one"),
    ("Sc1nccc2cn[nH]c12",  "1H-pyrazolo[5,4-c]pyridin-7(6H)-thione"),
    # 1H-pyrazolo[4,5-c]pyridine C3-OH/SH
    ("Oc1n[nH]c2ccncc12",  "1H-pyrazolo[4,5-c]pyridin-3(2H)-one"),
    ("Sc1n[nH]c2ccncc12",  "1H-pyrazolo[4,5-c]pyridin-3(2H)-thione"),
    # 1H-pyrazolo[4,5-c]pyridine C4-OH/SH
    ("Oc1nccc2[nH]ncc12",  "1H-pyrazolo[4,5-c]pyridin-4(3H)-one"),
    ("Sc1nccc2[nH]ncc12",  "1H-pyrazolo[4,5-c]pyridin-4(3H)-thione"),
    # 1H-pyrazolo[4,5-c]pyridine C6-OH/SH
    ("Oc1cc2[nH]ncc2cn1",  "1H-pyrazolo[4,5-c]pyridin-6(5H)-one"),
    ("Sc1cc2[nH]ncc2cn1",  "1H-pyrazolo[4,5-c]pyridin-6(5H)-thione"),
    # 1H-imidazo[4,5-c]pyridine C2-OH/SH
    ("Oc1nc2ccncc2[nH]1",  "1H-imidazo[4,5-c]pyridin-2(3H)-one"),
    ("Sc1nc2ccncc2[nH]1",  "1H-imidazo[4,5-c]pyridin-2(3H)-thione"),
    # 1H-imidazo[4,5-c]pyridine C4-OH/SH
    ("Oc1nccc2nc[nH]c12",  "1H-imidazo[4,5-c]pyridin-4(3H)-one"),
    ("Sc1nccc2nc[nH]c12",  "1H-imidazo[4,5-c]pyridin-4(3H)-thione"),
    # 1H-imidazo[4,5-c]pyridine C6-OH/SH
    ("Oc1cc2nc[nH]c2cn1",  "1H-imidazo[4,5-c]pyridin-6(5H)-one"),
    ("Sc1cc2nc[nH]c2cn1",  "1H-imidazo[4,5-c]pyridin-6(5H)-thione"),
    # Regressions: parent rings unchanged
    ("c1cncc2[nH]ncc12",   "1H-pyrazolo[5,4-c]pyridine"),
    ("c1cncc2cn[nH]c12",   "1H-pyrazolo[4,5-c]pyridine"),
    ("c1cncc2[nH]cnc12",   "1H-imidazo[4,5-c]pyridine"),
])
def test_phase793_pyrazolo_imidazo_pyridine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
