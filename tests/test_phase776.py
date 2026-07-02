"""Phase 776: fused 5+6 bicyclic α-ol/thiol → tautomers (IUPAC 2013).

- 1H/2H-indazole C3 (alpha to N2) → 1H-indazol-3(2H)-one/thione
- indolizine C3 and C5 (alpha to N4 junction) → 3(4H)-one/thione, 5(4H)-one/thione
- thieno[2,3-b]pyridine C6 (alpha to N1) → 6(1H)-one/thione
- furo[2,3-b]pyridine C6 (alpha to N1) → 6(1H)-one/thione
- 1H-pyrrolo[3,2-b]pyridine C5 (alpha to N4) → 5(4H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-indazole C3-OH/SH
    ("Oc1n[nH]c2ccccc12",       "1H-indazol-3(2H)-one"),
    ("Sc1n[nH]c2ccccc12",       "1H-indazol-3(2H)-thione"),
    # 2H-indazole C3-OH/SH (same preferred tautomer)
    ("Oc1[nH]nc2ccccc12",       "1H-indazol-3(2H)-one"),
    ("Sc1[nH]nc2ccccc12",       "1H-indazol-3(2H)-thione"),
    # indolizine C3 (alpha to N4 junction)
    ("Oc1ccc2ccccn12",          "indolizin-3(4H)-one"),
    ("Sc1ccc2ccccn12",          "indolizin-3(4H)-thione"),
    # indolizine C5 (alpha to N4 junction)
    ("Oc1cccc2cccn12",          "indolizin-5(4H)-one"),
    ("Sc1cccc2cccn12",          "indolizin-5(4H)-thione"),
    # thieno[2,3-b]pyridine C6 (alpha to N1)
    ("Oc1ccc2ccsc2n1",          "thieno[2,3-b]pyridin-6(1H)-one"),
    ("Sc1ccc2ccsc2n1",          "thieno[2,3-b]pyridin-6(1H)-thione"),
    # furo[2,3-b]pyridine C6 (alpha to N1)
    ("Oc1ccc2ccoc2n1",          "furo[2,3-b]pyridin-6(1H)-one"),
    ("Sc1ccc2ccoc2n1",          "furo[2,3-b]pyridin-6(1H)-thione"),
    # 1H-pyrrolo[3,2-b]pyridine C5 (alpha to N4)
    ("Oc1ccc2[nH]ccc2n1",       "1H-pyrrolo[3,2-b]pyridin-5(4H)-one"),
    ("Sc1ccc2[nH]ccc2n1",       "1H-pyrrolo[3,2-b]pyridin-5(4H)-thione"),
    # Regressions: parent rings unchanged
    ("c1ccc2[nH]ncc2c1",        "1H-indazole"),
    ("c1ccc2n[nH]cc2c1",        "2H-indazole"),
    ("c1ccn2cccc2c1",           "indolizine"),
    ("c1cnc2cc[nH]c2c1",        "1H-pyrrolo[3,2-b]pyridine"),
    # Regression: keto form input unchanged (Phase 417)
    ("O=c1[nH][nH]c2ccccc12",   "1H-indazol-3(2H)-one"),
])
def test_phase776_fused_5_6_bicyclic_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
