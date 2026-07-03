"""Phase 791: isothiazolo/isoxazolo[3,4-e]pyrazine, 1H-pyrazolo[4,5-b]pyridine α-ol/thiol → tautomers.

- isothiazolo[3,4-e]pyrazine C5 → 5(4H)-one/thione; C6 → 6(7H)-one/thione
- isoxazolo[3,4-e]pyrazine C5 → 5(4H)-one/thione; C6 → 6(7H)-one/thione
- 1H-pyrazolo[4,5-b]pyridine C3 → 3(2H)-one/thione; C5 → 5(4H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isothiazolo[3,4-e]pyrazine C5-OH/SH
    ("Oc1cnc2nscc2n1",    "isothiazolo[3,4-e]pyrazin-5(4H)-one"),
    ("Sc1cnc2nscc2n1",    "isothiazolo[3,4-e]pyrazin-5(4H)-thione"),
    # isothiazolo[3,4-e]pyrazine C6-OH/SH
    ("Oc1cnc2csnc2n1",    "isothiazolo[3,4-e]pyrazin-6(7H)-one"),
    ("Sc1cnc2csnc2n1",    "isothiazolo[3,4-e]pyrazin-6(7H)-thione"),
    # isoxazolo[3,4-e]pyrazine C5-OH/SH
    ("Oc1cnc2nocc2n1",    "isoxazolo[3,4-e]pyrazin-5(4H)-one"),
    ("Sc1cnc2nocc2n1",    "isoxazolo[3,4-e]pyrazin-5(4H)-thione"),
    # isoxazolo[3,4-e]pyrazine C6-OH/SH
    ("Oc1cnc2conc2n1",    "isoxazolo[3,4-e]pyrazin-6(7H)-one"),
    ("Sc1cnc2conc2n1",    "isoxazolo[3,4-e]pyrazin-6(7H)-thione"),
    # 1H-pyrazolo[4,5-b]pyridine C3-OH/SH
    ("Oc1n[nH]c2cccnc12", "1H-pyrazolo[4,5-b]pyridin-3(2H)-one"),
    ("Sc1n[nH]c2cccnc12", "1H-pyrazolo[4,5-b]pyridin-3(2H)-thione"),
    # 1H-pyrazolo[4,5-b]pyridine C5-OH/SH
    ("Oc1ccc2[nH]ncc2n1", "1H-pyrazolo[4,5-b]pyridin-5(4H)-one"),
    ("Sc1ccc2[nH]ncc2n1", "1H-pyrazolo[4,5-b]pyridin-5(4H)-thione"),
    # Regressions: parent rings unchanged
    ("c1cnc2csnc2n1",     "isothiazolo[3,4-e]pyrazine"),
    ("c1cnc2conc2n1",     "isoxazolo[3,4-e]pyrazine"),
    ("c1ccc2[nH]ncc2n1",  "1H-pyrazolo[4,5-b]pyridine"),
])
def test_phase791_fused_bicyclic_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
