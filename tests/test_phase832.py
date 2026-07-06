"""Phase 832: 1H-pyrazolo[3,4-c/d/e] α-ol/thiol → tautomers."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("[H]Oc1[nH]nc2nnccc12",         "pyrazolo[3,4-c]pyridazin-3(1H)-one"),
    ("[H]Sc1[nH]nc2nnccc12",         "pyrazolo[3,4-c]pyridazin-3(1H)-thione"),
    ("[H]Oc1cc2c[nH]nc2nn1",         "1H-pyrazolo[3,4-c]pyridazin-5(3H)-one"),
    ("[H]Sc1cc2c[nH]nc2nn1",         "1H-pyrazolo[3,4-c]pyridazin-5(3H)-thione"),
    ("[H]Oc1[nH]nc2cnncc12",         "pyrazolo[3,4-d]pyridazin-3(1H)-one"),
    ("[H]Sc1[nH]nc2cnncc12",         "pyrazolo[3,4-d]pyridazin-3(1H)-thione"),
    ("[H]Oc1nncc2n[nH]cc12",         "pyrazolo[3,4-d]pyridazin-4(1H)-one"),
    ("[H]Sc1nncc2n[nH]cc12",         "pyrazolo[3,4-d]pyridazin-4(1H)-thione"),
    ("[H]Oc1nncc2c[nH]nc12",         "1H-pyrazolo[3,4-d]pyridazin-7(3H)-one"),
    ("[H]Sc1nncc2c[nH]nc12",         "1H-pyrazolo[3,4-d]pyridazin-7(3H)-thione"),
    ("[H]Oc1nnnc2n[nH]cc12",         "1H-pyrazolo[3,4-d][1,2,3]triazin-4(4H)-one"),
    ("[H]Sc1nnnc2n[nH]cc12",         "1H-pyrazolo[3,4-d][1,2,3]triazin-4(4H)-thione"),
    ("[H]Oc1[nH]nc2nnncc12",         "pyrazolo[3,4-d][1,2,3]triazin-5(1H)-one"),
    ("[H]Sc1[nH]nc2nnncc12",         "pyrazolo[3,4-d][1,2,3]triazin-5(1H)-thione"),
    ("[H]Oc1n[nH]c2cnnc-2n1",        "pyrazolo[3,4-e][1,2,4]triazin-3(1H)-one"),
    ("[H]Sc1n[nH]c2cnnc-2n1",        "pyrazolo[3,4-e][1,2,4]triazin-3(1H)-thione"),
    ("[H]Oc1nnc2ncn[nH]c1-2",        "pyrazolo[3,4-e][1,2,4]triazin-7(1H)-one"),
    ("[H]Sc1nnc2ncn[nH]c1-2",        "pyrazolo[3,4-e][1,2,4]triazin-7(1H)-thione"),
    ("[H]Oc1[nH]nc2nccnc12",         "pyrazolo[3,4-e]pyrazin-3(1H)-one"),
    ("[H]Sc1[nH]nc2nccnc12",         "pyrazolo[3,4-e]pyrazin-3(1H)-thione"),
    ("[H]Oc1cnc2n[nH]cc2n1",         "1H-pyrazolo[3,4-e]pyrazin-5(3H)-one"),
    ("[H]Sc1cnc2n[nH]cc2n1",         "1H-pyrazolo[3,4-e]pyrazin-5(3H)-thione"),
    ("[H]Oc1cnc2c[nH]nc2n1",         "pyrazolo[3,4-e]pyrazin-6(1H)-one"),
    ("[H]Sc1cnc2c[nH]nc2n1",         "pyrazolo[3,4-e]pyrazin-6(1H)-thione"),
    ("c1cc2c[nH]nc2nn1",              "1H-pyrazolo[3,4-c]pyridazine"),
    ("c1nncc2n[nH]cc12",             "1H-pyrazolo[3,4-d]pyridazine"),
    ("c1nnnc2n[nH]cc12",             "1H-pyrazolo[3,4-d][1,2,3]triazine"),
    ("c1n[nH]c2cnnc-2n1",            "1H-pyrazolo[3,4-e][1,2,4]triazine"),
    ("c1cnc2n[nH]cc2n1",             "1H-pyrazolo[3,4-e]pyrazine"),
])
def test_phase832(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
