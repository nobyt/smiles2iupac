"""Phase 829: 1H-pyrrolo[2,3-c/d/e] α-ol/thiol → tautomers."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("[H]Oc1cc2cc[nH]c2nn1",      "4H-pyrrolo[2,3-c]pyridazin-3(7H)-one"),
    ("[H]Sc1cc2cc[nH]c2nn1",      "4H-pyrrolo[2,3-c]pyridazin-3(7H)-thione"),
    ("[H]Oc1cc2ccnnc2[nH]1",      "5H-pyrrolo[2,3-c]pyridazin-6(7H)-one"),
    ("[H]Sc1cc2ccnnc2[nH]1",      "5H-pyrrolo[2,3-c]pyridazin-6(7H)-thione"),
    ("[H]Oc1cc2cnncc2[nH]1",      "pyrrolo[2,3-d]pyridazin-2(1H)-one"),
    ("[H]Sc1cc2cnncc2[nH]1",      "pyrrolo[2,3-d]pyridazin-2(1H)-thione"),
    ("[H]Oc1nncc2[nH]ccc12",      "pyrrolo[2,3-d]pyridazin-4(1H)-one"),
    ("[H]Sc1nncc2[nH]ccc12",      "pyrrolo[2,3-d]pyridazin-4(1H)-thione"),
    ("[H]Oc1nncc2cc[nH]c12",      "1H-pyrrolo[2,3-d]pyridazin-7(3H)-one"),
    ("[H]Sc1nncc2cc[nH]c12",      "1H-pyrrolo[2,3-d]pyridazin-7(3H)-thione"),
    ("[H]Oc1ncc2cc[nH]c2n1",      "1H-pyrrolo[2,3-d]pyrimidin-2(7H)-one"),
    ("[H]Sc1ncc2cc[nH]c2n1",      "1H-pyrrolo[2,3-d]pyrimidin-2(7H)-thione"),
    ("[H]Oc1ncnc2[nH]ccc12",      "1H-pyrrolo[2,3-d]pyrimidin-4(7H)-one"),
    ("[H]Sc1ncnc2[nH]ccc12",      "1H-pyrrolo[2,3-d]pyrimidin-4(7H)-thione"),
    ("[H]Oc1cc2cncnc2[nH]1",      "5H-pyrrolo[2,3-d]pyrimidin-6(7H)-one"),
    ("[H]Sc1cc2cncnc2[nH]1",      "5H-pyrrolo[2,3-d]pyrimidin-6(7H)-thione"),
    ("[H]Oc1nnnc2[nH]ccc12",      "1H-pyrrolo[2,3-d][1,2,3]triazin-4(7H)-one"),
    ("[H]Sc1nnnc2[nH]ccc12",      "1H-pyrrolo[2,3-d][1,2,3]triazin-4(7H)-thione"),
    ("[H]Oc1cc2cnnnc2[nH]1",      "5H-pyrrolo[2,3-d][1,2,3]triazin-6(7H)-one"),
    ("[H]Sc1cc2cnnnc2[nH]1",      "5H-pyrrolo[2,3-d][1,2,3]triazin-6(7H)-thione"),
    ("[H]Oc1n[nH]c2ccnc-2n1",    "pyrrolo[2,3-e][1,2,4]triazin-3(1H)-one"),
    ("[H]Sc1n[nH]c2ccnc-2n1",    "pyrrolo[2,3-e][1,2,4]triazin-3(1H)-thione"),
    ("[H]Oc1cc2[nH]ncnc-2n1",    "pyrrolo[2,3-e][1,2,4]triazin-6(1H)-one"),
    ("[H]Sc1cc2[nH]ncnc-2n1",    "pyrrolo[2,3-e][1,2,4]triazin-6(1H)-thione"),
    ("[H]Oc1cnc2nccc-2[nH]1",    "pyrrolo[2,3-e]pyrazin-2(1H)-one"),
    ("[H]Sc1cnc2nccc-2[nH]1",    "pyrrolo[2,3-e]pyrazin-2(1H)-thione"),
    ("[H]Oc1c[nH]c2ccnc-2n1",    "pyrrolo[2,3-e]pyrazin-3(1H)-one"),
    ("[H]Sc1c[nH]c2ccnc-2n1",    "pyrrolo[2,3-e]pyrazin-3(1H)-thione"),
    ("[H]Oc1cc2[nH]ccnc-2n1",    "pyrrolo[2,3-e]pyrazin-6(1H)-one"),
    ("[H]Sc1cc2[nH]ccnc-2n1",    "pyrrolo[2,3-e]pyrazin-6(1H)-thione"),
    ("c1cc2cc[nH]c2nn1",          "7H-pyrrolo[2,3-c]pyridazine"),
    ("c1cc2cnncc2[nH]1",          "1H-pyrrolo[2,3-d]pyridazine"),
    ("c1ncc2cc[nH]c2n1",          "7H-pyrrolo[2,3-d]pyrimidine"),
    ("c1cc2cnnnc2[nH]1",          "7H-pyrrolo[2,3-d][1,2,3]triazine"),
    ("c1cc2[nH]ncnc-2n1",         "1H-pyrrolo[2,3-e][1,2,4]triazine"),
    ("c1c[nH]c2ccnc-2n1",         "1H-pyrrolo[2,3-e]pyrazine"),
])
def test_phase829(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
