"""Phase 790: thieno/furo[3,2-e][1,2,4]triazine, isothiazolo/isoxazolo[5,4-c]pyridazine, pyrazolo[1,5-b]pyridazine, pyrazolo[1,5-a]pyridine C7 α-ol/thiol → tautomers.

- thieno[3,2-e][1,2,4]triazine C3 → 3(2H)-one/thione
- furo[3,2-e][1,2,4]triazine C3 → 3(2H)-one/thione
- isothiazolo[5,4-c]pyridazine C3 → 3(4H)-one/thione; C5 → 5(4H)-one/thione
- isoxazolo[5,4-c]pyridazine C3 → 3(4H)-one/thione; C5 → 5(4H)-one/thione
- pyrazolo[1,5-b]pyridazine C2 → 2(1H)-one/thione; C6 → 6(5H)-one/thione
- pyrazolo[1,5-a]pyridine C7 → 7(4H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[3,2-e][1,2,4]triazine C3-OH/SH
    ("Oc1nnc2sccc2n1",    "thieno[3,2-e][1,2,4]triazin-3(2H)-one"),
    ("Sc1nnc2sccc2n1",    "thieno[3,2-e][1,2,4]triazin-3(2H)-thione"),
    # furo[3,2-e][1,2,4]triazine C3-OH/SH
    ("Oc1nnc2occc2n1",    "furo[3,2-e][1,2,4]triazin-3(2H)-one"),
    ("Sc1nnc2occc2n1",    "furo[3,2-e][1,2,4]triazin-3(2H)-thione"),
    # isothiazolo[5,4-c]pyridazine C3-OH/SH
    ("Oc1nsc2nnccc12",    "isothiazolo[5,4-c]pyridazin-3(4H)-one"),
    ("Sc1nsc2nnccc12",    "isothiazolo[5,4-c]pyridazin-3(4H)-thione"),
    # isothiazolo[5,4-c]pyridazine C5-OH/SH
    ("Oc1cc2cnsc2nn1",    "isothiazolo[5,4-c]pyridazin-5(4H)-one"),
    ("Sc1cc2cnsc2nn1",    "isothiazolo[5,4-c]pyridazin-5(4H)-thione"),
    # isoxazolo[5,4-c]pyridazine C3-OH/SH
    ("Oc1noc2nnccc12",    "isoxazolo[5,4-c]pyridazin-3(4H)-one"),
    ("Sc1noc2nnccc12",    "isoxazolo[5,4-c]pyridazin-3(4H)-thione"),
    # isoxazolo[5,4-c]pyridazine C5-OH/SH
    ("Oc1cc2cnoc2nn1",    "isoxazolo[5,4-c]pyridazin-5(4H)-one"),
    ("Sc1cc2cnoc2nn1",    "isoxazolo[5,4-c]pyridazin-5(4H)-thione"),
    # pyrazolo[1,5-b]pyridazine C2-OH/SH
    ("Oc1cc2cccnn2n1",    "pyrazolo[1,5-b]pyridazin-2(1H)-one"),
    ("Sc1cc2cccnn2n1",    "pyrazolo[1,5-b]pyridazin-2(1H)-thione"),
    # pyrazolo[1,5-b]pyridazine C6-OH/SH
    ("Oc1ccc2ccnn2n1",    "pyrazolo[1,5-b]pyridazin-6(5H)-one"),
    ("Sc1ccc2ccnn2n1",    "pyrazolo[1,5-b]pyridazin-6(5H)-thione"),
    # pyrazolo[1,5-a]pyridine C7-OH/SH
    ("Oc1cccc2ccnn12",    "pyrazolo[1,5-a]pyridin-7(4H)-one"),
    ("Sc1cccc2ccnn12",    "pyrazolo[1,5-a]pyridin-7(4H)-thione"),
    # Regressions: parent rings unchanged
    ("c1csc2nncnc12",     "thieno[3,2-e][1,2,4]triazine"),
    ("c1coc2nncnc12",     "furo[3,2-e][1,2,4]triazine"),
    ("c1nsc2nnccc12",     "isothiazolo[5,4-c]pyridazine"),
    ("c1noc2nnccc12",     "isoxazolo[5,4-c]pyridazine"),
    ("c1cnn2ncccc12",     "pyrazolo[1,5-b]pyridazine"),
])
def test_phase790_fused_bicyclic_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
