"""Phase 785: thieno/furo[2,3-e]pyrazine, isothiazolo/isoxazolo[5,4-b]pyridine, thieno[3,2-d]pyrimidine α-ol/thiol → tautomers.

- thieno[2,3-e]pyrazine C2 → 2(1H)-one/thione; C3 → 3(4H)-one/thione
- furo[2,3-e]pyrazine C2 → 2(1H)-one/thione; C3 → 3(4H)-one/thione
- isothiazolo[5,4-b]pyridine C3 → 3(4H)-one/thione; C6 → 6(5H)-one/thione
- isoxazolo[5,4-b]pyridine C3 → 3(4H)-one/thione; C6 → 6(5H)-one/thione
- thieno[3,2-d]pyrimidine C2 → 2(1H)-one/thione; C4 → 4(3H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[2,3-e]pyrazine C2-OH/SH
    ("Oc1cnc2sccc2n1",    "thieno[2,3-e]pyrazin-2(1H)-one"),
    ("Sc1cnc2sccc2n1",    "thieno[2,3-e]pyrazin-2(1H)-thione"),
    # thieno[2,3-e]pyrazine C3-OH/SH
    ("Oc1cnc2ccsc2n1",    "thieno[2,3-e]pyrazin-3(4H)-one"),
    ("Sc1cnc2ccsc2n1",    "thieno[2,3-e]pyrazin-3(4H)-thione"),
    # furo[2,3-e]pyrazine C2-OH/SH
    ("Oc1cnc2occc2n1",    "furo[2,3-e]pyrazin-2(1H)-one"),
    ("Sc1cnc2occc2n1",    "furo[2,3-e]pyrazin-2(1H)-thione"),
    # furo[2,3-e]pyrazine C3-OH/SH
    ("Oc1cnc2ccoc2n1",    "furo[2,3-e]pyrazin-3(4H)-one"),
    ("Sc1cnc2ccoc2n1",    "furo[2,3-e]pyrazin-3(4H)-thione"),
    # isothiazolo[5,4-b]pyridine C3-OH/SH
    ("Oc1nsc2ncccc12",    "isothiazolo[5,4-b]pyridin-3(4H)-one"),
    ("Sc1nsc2ncccc12",    "isothiazolo[5,4-b]pyridin-3(4H)-thione"),
    # isothiazolo[5,4-b]pyridine C6-OH/SH
    ("Oc1ccc2cnsc2n1",    "isothiazolo[5,4-b]pyridin-6(5H)-one"),
    ("Sc1ccc2cnsc2n1",    "isothiazolo[5,4-b]pyridin-6(5H)-thione"),
    # isoxazolo[5,4-b]pyridine C3-OH/SH
    ("Oc1noc2ncccc12",    "isoxazolo[5,4-b]pyridin-3(4H)-one"),
    ("Sc1noc2ncccc12",    "isoxazolo[5,4-b]pyridin-3(4H)-thione"),
    # isoxazolo[5,4-b]pyridine C6-OH/SH
    ("Oc1ccc2cnoc2n1",    "isoxazolo[5,4-b]pyridin-6(5H)-one"),
    ("Sc1ccc2cnoc2n1",    "isoxazolo[5,4-b]pyridin-6(5H)-thione"),
    # thieno[3,2-d]pyrimidine C2-OH/SH
    ("Oc1ncc2sccc2n1",    "thieno[3,2-d]pyrimidin-2(1H)-one"),
    ("Sc1ncc2sccc2n1",    "thieno[3,2-d]pyrimidin-2(1H)-thione"),
    # thieno[3,2-d]pyrimidine C4-OH/SH
    ("Oc1ncnc2ccsc12",    "thieno[3,2-d]pyrimidin-4(3H)-one"),
    ("Sc1ncnc2ccsc12",    "thieno[3,2-d]pyrimidin-4(3H)-thione"),
    # Regressions: parent rings unchanged
    ("c1csc2nccnc12",     "thieno[2,3-e]pyrazine"),
    ("c1coc2nccnc12",     "furo[2,3-e]pyrazine"),
    ("c1nsc2ncccc12",     "isothiazolo[5,4-b]pyridine"),
    ("c1noc2ncccc12",     "isoxazolo[5,4-b]pyridine"),
    ("c1csc2cncnc12",     "thieno[3,2-d]pyrimidine"),
])
def test_phase785_thienofuro_pyrazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
