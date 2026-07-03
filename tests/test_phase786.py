"""Phase 786: furo[3,2-d]pyrimidine, isothiazolo/isoxazolo[4,5-d]pyrimidine, 1H-indole α-ol/thiol → tautomers.

- furo[3,2-d]pyrimidine C2 → 2(1H)-one/thione; C4 → 4(3H)-one/thione
- isothiazolo[4,5-d]pyrimidine C3 → 3(2H)-one/thione; C5 → 5(4H)-one/thione; C7 → 7(6H)-one/thione
- isoxazolo[4,5-d]pyrimidine C3 → 3(2H)-one/thione; C5 → 5(4H)-one/thione; C7 → 7(6H)-one/thione
- 1H-indole C2 → 1H-indol-2(3H)-one/thione (oxindole/thiooxindole)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # furo[3,2-d]pyrimidine C2-OH/SH
    ("Oc1ncc2occc2n1",    "furo[3,2-d]pyrimidin-2(1H)-one"),
    ("Sc1ncc2occc2n1",    "furo[3,2-d]pyrimidin-2(1H)-thione"),
    # furo[3,2-d]pyrimidine C4-OH/SH
    ("Oc1ncnc2ccoc12",    "furo[3,2-d]pyrimidin-4(3H)-one"),
    ("Sc1ncnc2ccoc12",    "furo[3,2-d]pyrimidin-4(3H)-thione"),
    # isothiazolo[4,5-d]pyrimidine C3-OH/SH
    ("Oc1nsc2cncnc12",    "isothiazolo[4,5-d]pyrimidin-3(2H)-one"),
    ("Sc1nsc2cncnc12",    "isothiazolo[4,5-d]pyrimidin-3(2H)-thione"),
    # isothiazolo[4,5-d]pyrimidine C5-OH/SH
    ("Oc1ncc2sncc2n1",    "isothiazolo[4,5-d]pyrimidin-5(4H)-one"),
    ("Sc1ncc2sncc2n1",    "isothiazolo[4,5-d]pyrimidin-5(4H)-thione"),
    # isothiazolo[4,5-d]pyrimidine C7-OH/SH
    ("Oc1ncnc2cnsc12",    "isothiazolo[4,5-d]pyrimidin-7(6H)-one"),
    ("Sc1ncnc2cnsc12",    "isothiazolo[4,5-d]pyrimidin-7(6H)-thione"),
    # isoxazolo[4,5-d]pyrimidine C3-OH/SH
    ("Oc1noc2cncnc12",    "isoxazolo[4,5-d]pyrimidin-3(2H)-one"),
    ("Sc1noc2cncnc12",    "isoxazolo[4,5-d]pyrimidin-3(2H)-thione"),
    # isoxazolo[4,5-d]pyrimidine C5-OH/SH
    ("Oc1ncc2oncc2n1",    "isoxazolo[4,5-d]pyrimidin-5(4H)-one"),
    ("Sc1ncc2oncc2n1",    "isoxazolo[4,5-d]pyrimidin-5(4H)-thione"),
    # isoxazolo[4,5-d]pyrimidine C7-OH/SH
    ("Oc1ncnc2cnoc12",    "isoxazolo[4,5-d]pyrimidin-7(6H)-one"),
    ("Sc1ncnc2cnoc12",    "isoxazolo[4,5-d]pyrimidin-7(6H)-thione"),
    # 1H-indole C2-OH/SH
    ("Oc1cc2ccccc2[nH]1", "1H-indol-2(3H)-one"),
    ("Sc1cc2ccccc2[nH]1", "1H-indol-2(3H)-thione"),
    # Regressions: parent rings unchanged
    ("c1coc2cncnc12",     "furo[3,2-d]pyrimidine"),
    ("c1nsc2cncnc12",     "isothiazolo[4,5-d]pyrimidine"),
    ("c1noc2cncnc12",     "isoxazolo[4,5-d]pyrimidine"),
    ("c1ccc2[nH]ccc2c1",  "1H-indole"),
])
def test_phase786_furo_indole_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
