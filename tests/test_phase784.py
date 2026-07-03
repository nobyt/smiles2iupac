"""Phase 784: thieno/furo[2,3-d]pyrimidine, isothiazolo/isoxazolo[5,4-d]pyrimidine, thieno[3,2-b]pyridine α-ol/thiol → tautomers.

- thieno[2,3-d]pyrimidine C2 → 2(1H)-one/thione; C4 → 4(3H)-one/thione
- furo[2,3-d]pyrimidine C2 → 2(1H)-one/thione; C4 → 4(3H)-one/thione
- isothiazolo[5,4-d]pyrimidine C3 → 3(4H)-one/thione; C4 → 4(3H)-one/thione; C6 → 6(5H)-one/thione
- isoxazolo[5,4-d]pyrimidine C3 → 3(4H)-one/thione; C4 → 4(3H)-one/thione; C6 → 6(5H)-one/thione
- thieno[3,2-b]pyridine C5 → 5(4H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[2,3-d]pyrimidine C2-OH/SH
    ("Oc1ncc2ccsc2n1",    "thieno[2,3-d]pyrimidin-2(1H)-one"),
    ("Sc1ncc2ccsc2n1",    "thieno[2,3-d]pyrimidin-2(1H)-thione"),
    # thieno[2,3-d]pyrimidine C4-OH/SH
    ("Oc1ncnc2sccc12",    "thieno[2,3-d]pyrimidin-4(3H)-one"),
    ("Sc1ncnc2sccc12",    "thieno[2,3-d]pyrimidin-4(3H)-thione"),
    # furo[2,3-d]pyrimidine C2-OH/SH
    ("Oc1ncc2ccoc2n1",    "furo[2,3-d]pyrimidin-2(1H)-one"),
    ("Sc1ncc2ccoc2n1",    "furo[2,3-d]pyrimidin-2(1H)-thione"),
    # furo[2,3-d]pyrimidine C4-OH/SH
    ("Oc1ncnc2occc12",    "furo[2,3-d]pyrimidin-4(3H)-one"),
    ("Sc1ncnc2occc12",    "furo[2,3-d]pyrimidin-4(3H)-thione"),
    # isothiazolo[5,4-d]pyrimidine C3-OH/SH
    ("Oc1nsc2ncncc12",    "isothiazolo[5,4-d]pyrimidin-3(4H)-one"),
    ("Sc1nsc2ncncc12",    "isothiazolo[5,4-d]pyrimidin-3(4H)-thione"),
    # isothiazolo[5,4-d]pyrimidine C4-OH/SH
    ("Oc1ncnc2sncc12",    "isothiazolo[5,4-d]pyrimidin-4(3H)-one"),
    ("Sc1ncnc2sncc12",    "isothiazolo[5,4-d]pyrimidin-4(3H)-thione"),
    # isothiazolo[5,4-d]pyrimidine C6-OH/SH
    ("Oc1ncc2cnsc2n1",    "isothiazolo[5,4-d]pyrimidin-6(5H)-one"),
    ("Sc1ncc2cnsc2n1",    "isothiazolo[5,4-d]pyrimidin-6(5H)-thione"),
    # isoxazolo[5,4-d]pyrimidine C3-OH/SH
    ("Oc1noc2ncncc12",    "isoxazolo[5,4-d]pyrimidin-3(4H)-one"),
    ("Sc1noc2ncncc12",    "isoxazolo[5,4-d]pyrimidin-3(4H)-thione"),
    # isoxazolo[5,4-d]pyrimidine C4-OH/SH
    ("Oc1ncnc2oncc12",    "isoxazolo[5,4-d]pyrimidin-4(3H)-one"),
    ("Sc1ncnc2oncc12",    "isoxazolo[5,4-d]pyrimidin-4(3H)-thione"),
    # isoxazolo[5,4-d]pyrimidine C6-OH/SH
    ("Oc1ncc2cnoc2n1",    "isoxazolo[5,4-d]pyrimidin-6(5H)-one"),
    ("Sc1ncc2cnoc2n1",    "isoxazolo[5,4-d]pyrimidin-6(5H)-thione"),
    # thieno[3,2-b]pyridine C5-OH/SH
    ("Oc1ccc2sccc2n1",    "thieno[3,2-b]pyridin-5(4H)-one"),
    ("Sc1ccc2sccc2n1",    "thieno[3,2-b]pyridin-5(4H)-thione"),
    # Regressions: parent rings unchanged
    ("c1csc2ncncc12",     "thieno[2,3-d]pyrimidine"),
    ("c1coc2ncncc12",     "furo[2,3-d]pyrimidine"),
    ("c1nsc2ncncc12",     "isothiazolo[5,4-d]pyrimidine"),
    ("c1noc2ncncc12",     "isoxazolo[5,4-d]pyrimidine"),
    ("c1csc2cccnc12",     "thieno[3,2-b]pyridine"),
])
def test_phase784_thienofuro_pyrimidine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
