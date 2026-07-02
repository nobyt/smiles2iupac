"""Phase 783: pyrido[2,3-c]pyridazine, isothiazolo[4,3-b/4,5-b]pyridine, 1H-pyrrolo[3,2-b/2,3-b]pyridine α-ol/thiol → tautomers.

- pyrido[2,3-c]pyridazine C3-ol/thiol → 3(2H)-one/thione; C7-ol/thiol → 7(6H)-one/thione
- isothiazolo[4,3-b]pyridine C5-ol/thiol → 5(4H)-one/thione
- isothiazolo[4,5-b]pyridine C5-ol/thiol → 5(4H)-one/thione; C3-ol/thiol → 3(4H)-one/thione
- 1H-pyrrolo[3,2-b]pyridine C2-ol/thiol → 2(3H)-one/thione
- 1H-pyrrolo[2,3-b]pyridine C2-ol/thiol → 2(3H)-one/thione; C6-ol/thiol → 6(5H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrido[2,3-c]pyridazine C3-OH/SH
    ("Oc1cc2cccnc2nn1",    "pyrido[2,3-c]pyridazin-3(2H)-one"),
    ("Sc1cc2cccnc2nn1",    "pyrido[2,3-c]pyridazin-3(2H)-thione"),
    # pyrido[2,3-c]pyridazine C7-OH/SH
    ("Oc1ccc2ccnnc2n1",    "pyrido[2,3-c]pyridazin-7(6H)-one"),
    ("Sc1ccc2ccnnc2n1",    "pyrido[2,3-c]pyridazin-7(6H)-thione"),
    # isothiazolo[4,3-b]pyridine C5-OH/SH
    ("Oc1ccc2nscc2n1",     "isothiazolo[4,3-b]pyridin-5(4H)-one"),
    ("Sc1ccc2nscc2n1",     "isothiazolo[4,3-b]pyridin-5(4H)-thione"),
    # isothiazolo[4,5-b]pyridine C5-OH/SH
    ("Oc1ccc2sncc2n1",     "isothiazolo[4,5-b]pyridin-5(4H)-one"),
    ("Sc1ccc2sncc2n1",     "isothiazolo[4,5-b]pyridin-5(4H)-thione"),
    # isothiazolo[4,5-b]pyridine C3-OH/SH
    ("Oc1nsc2cccnc12",     "isothiazolo[4,5-b]pyridin-3(4H)-one"),
    ("Sc1nsc2cccnc12",     "isothiazolo[4,5-b]pyridin-3(4H)-thione"),
    # 1H-pyrrolo[3,2-b]pyridine C2-OH/SH
    ("Oc1cc2ncccc2[nH]1",  "1H-pyrrolo[3,2-b]pyridin-2(3H)-one"),
    ("Sc1cc2ncccc2[nH]1",  "1H-pyrrolo[3,2-b]pyridin-2(3H)-thione"),
    # 1H-pyrrolo[2,3-b]pyridine C2-OH/SH
    ("Oc1cc2cccnc2[nH]1",  "1H-pyrrolo[2,3-b]pyridin-2(3H)-one"),
    ("Sc1cc2cccnc2[nH]1",  "1H-pyrrolo[2,3-b]pyridin-2(3H)-thione"),
    # 1H-pyrrolo[2,3-b]pyridine C6-OH/SH
    ("Oc1ccc2cc[nH]c2n1",  "1H-pyrrolo[2,3-b]pyridin-6(5H)-one"),
    ("Sc1ccc2cc[nH]c2n1",  "1H-pyrrolo[2,3-b]pyridin-6(5H)-thione"),
    # Regressions: parent rings unchanged
    ("c1cnc2nnccc2c1",     "pyrido[2,3-c]pyridazine"),
    ("c1cnc2csnc2c1",      "isothiazolo[4,3-b]pyridine"),
    ("c1cnc2cnsc2c1",      "isothiazolo[4,5-b]pyridine"),
    ("c1cnc2cc[nH]c2c1",   "1H-pyrrolo[3,2-b]pyridine"),
    ("c1cnc2[nH]ccc2c1",   "1H-pyrrolo[2,3-b]pyridine"),
])
def test_phase783_fused_bicyclic_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
