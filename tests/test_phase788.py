"""Phase 788: 9H-purine, 3H-imidazo[4,5-b]pyridine, 1H-pyrazolo[5,4-d]pyrimidine α-ol/thiol → tautomers.

- 9H-purine C8 → 9H-purin-8(7H)-one/thione
- 3H-imidazo[4,5-b]pyridine C2 → 1H-imidazo[4,5-b]pyridin-2(3H)-one/thione (preferred PIN)
- 3H-imidazo[4,5-b]pyridine C5 → 1H-imidazo[4,5-b]pyridin-5(4H)-one/thione (preferred PIN)
- 1H-pyrazolo[5,4-d]pyrimidine C3 → 3(2H)-one/thione; C4 → 4(3H)-one/thione; C6 → 6(5H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 9H-purine C8-OH/SH
    ("Oc1nc2cncnc2[nH]1",   "9H-purin-8(7H)-one"),
    ("Sc1nc2cncnc2[nH]1",   "9H-purin-8(7H)-thione"),
    # 3H-imidazo[4,5-b]pyridine C2-OH/SH → preferred PIN 1H-form
    ("Oc1nc2ncccc2[nH]1",   "1H-imidazo[4,5-b]pyridin-2(3H)-one"),
    ("Sc1nc2ncccc2[nH]1",   "1H-imidazo[4,5-b]pyridin-2(3H)-thione"),
    # 3H-imidazo[4,5-b]pyridine C5-OH/SH → preferred PIN 1H-form
    ("Oc1ccc2[nH]cnc2n1",   "1H-imidazo[4,5-b]pyridin-5(4H)-one"),
    ("Sc1ccc2[nH]cnc2n1",   "1H-imidazo[4,5-b]pyridin-5(4H)-thione"),
    # 1H-pyrazolo[5,4-d]pyrimidine C3-OH/SH
    ("Oc1n[nH]c2ncncc12",   "1H-pyrazolo[5,4-d]pyrimidin-3(2H)-one"),
    ("Sc1n[nH]c2ncncc12",   "1H-pyrazolo[5,4-d]pyrimidin-3(2H)-thione"),
    # 1H-pyrazolo[5,4-d]pyrimidine C4-OH/SH
    ("Oc1ncnc2[nH]ncc12",   "1H-pyrazolo[5,4-d]pyrimidin-4(3H)-one"),
    ("Sc1ncnc2[nH]ncc12",   "1H-pyrazolo[5,4-d]pyrimidin-4(3H)-thione"),
    # 1H-pyrazolo[5,4-d]pyrimidine C6-OH/SH
    ("Oc1ncc2cn[nH]c2n1",   "1H-pyrazolo[5,4-d]pyrimidin-6(5H)-one"),
    ("Sc1ncc2cn[nH]c2n1",   "1H-pyrazolo[5,4-d]pyrimidin-6(5H)-thione"),
    # Regressions: parent rings unchanged
    ("c1nc2cncnc2[nH]1",    "9H-purine"),
    ("c1ccc2[nH]cnc2n1",    "3H-imidazo[4,5-b]pyridine"),
    ("c1n[nH]c2ncncc12",    "1H-pyrazolo[5,4-d]pyrimidine"),
])
def test_phase788_purine_imidazo_pyrazolo_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
