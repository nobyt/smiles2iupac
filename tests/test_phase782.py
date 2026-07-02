"""Phase 782: 1H/2H-[1,2,3]triazolo[4,5-b]pyridine and isothiazolo[3,4-b]pyridine α-ol/thiol → tautomers.

- 1H/2H-[1,2,3]triazolo[4,5-b]pyridine C5-ol/thiol → 1H-[1,2,3]triazolo[4,5-b]pyridin-5(4H)-one/thione
  (both tautomeric H-forms of the parent give the same preferred PIN)
- isothiazolo[3,4-b]pyridine C6-ol/thiol → isothiazolo[3,4-b]pyridin-6(5H)-one/thione
  (H on N5 because N1 is in fact S at position 1)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-[1,2,3]triazolo[4,5-b]pyridine C5-OH/SH
    ("Oc1ccc2nn[nH]c2n1",   "1H-[1,2,3]triazolo[4,5-b]pyridin-5(4H)-one"),
    ("Sc1ccc2nn[nH]c2n1",   "1H-[1,2,3]triazolo[4,5-b]pyridin-5(4H)-thione"),
    # 2H form gives same preferred tautomer
    ("Oc1ccc2n[nH]nc2n1",   "1H-[1,2,3]triazolo[4,5-b]pyridin-5(4H)-one"),
    ("Sc1ccc2n[nH]nc2n1",   "1H-[1,2,3]triazolo[4,5-b]pyridin-5(4H)-thione"),
    # isothiazolo[3,4-b]pyridine C6-OH/SH
    ("Oc1ccc2csnc2n1",       "isothiazolo[3,4-b]pyridin-6(5H)-one"),
    ("Sc1ccc2csnc2n1",       "isothiazolo[3,4-b]pyridin-6(5H)-thione"),
    # Regressions: parent rings unchanged
    ("c1cnc2[nH]nnc2c1",    "1H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("c1cnc2n[nH]nc2c1",    "2H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("c1cnc2nscc2c1",        "isothiazolo[3,4-b]pyridine"),
])
def test_phase782_triazolo_isothiazolo_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
