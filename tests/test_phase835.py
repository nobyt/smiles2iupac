"""Phase 835: 1H-imidazo[4,5-c/d/e] + 1H-imidazo[5,4-d/e] α-ol/thiol → tautomers."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("[H]Oc1cc2[nH]cnc2nn1",         "imidazo[4,5-c]pyridazin-3(1H)-one"),
    ("[H]Sc1cc2[nH]cnc2nn1",         "imidazo[4,5-c]pyridazin-3(1H)-thione"),
    ("[H]Oc1nc2nnccc2[nH]1",         "1H-imidazo[4,5-c]pyridazin-6(5H)-one"),
    ("[H]Sc1nc2nnccc2[nH]1",         "1H-imidazo[4,5-c]pyridazin-6(5H)-thione"),
    ("[H]Oc1nc2cnncc2[nH]1",         "1H-imidazo[4,5-d]pyridazin-2(1H)-one"),
    ("[H]Sc1nc2cnncc2[nH]1",         "1H-imidazo[4,5-d]pyridazin-2(1H)-thione"),
    ("[H]Oc1nncc2[nH]cnc12",         "1H-imidazo[4,5-d]pyridazin-4(1H)-one"),
    ("[H]Sc1nncc2[nH]cnc12",         "1H-imidazo[4,5-d]pyridazin-4(1H)-thione"),
    ("[H]Oc1nnnc2nc[nH]c12",         "1H-imidazo[4,5-d][1,2,3]triazin-4(5H)-one"),
    ("[H]Sc1nnnc2nc[nH]c12",         "1H-imidazo[4,5-d][1,2,3]triazin-4(5H)-thione"),
    ("[H]Oc1nc2nnncc2[nH]1",         "1H-imidazo[4,5-d][1,2,3]triazin-6(5H)-one"),
    ("[H]Sc1nc2nnncc2[nH]1",         "1H-imidazo[4,5-d][1,2,3]triazin-6(5H)-thione"),
    ("[H]Oc1nnc2[nH]cnc2n1",         "3H-imidazo[4,5-e][1,2,4]triazin-3(7H)-one"),
    ("[H]Sc1nnc2[nH]cnc2n1",         "3H-imidazo[4,5-e][1,2,4]triazin-3(7H)-thione"),
    ("[H]Oc1nc2ncnnc2[nH]1",         "4H-imidazo[4,5-e][1,2,4]triazin-6(7H)-one"),
    ("[H]Sc1nc2ncnnc2[nH]1",         "4H-imidazo[4,5-e][1,2,4]triazin-6(7H)-thione"),
    ("[H]Oc1nnnc2[nH]cnc12",         "1H-imidazo[4,5-d][1,2,3]triazin-4(5H)-one"),
    ("[H]Sc1nnnc2[nH]cnc12",         "1H-imidazo[4,5-d][1,2,3]triazin-4(5H)-thione"),
    ("[H]Oc1nc2cnnnc2[nH]1",         "1H-imidazo[4,5-d][1,2,3]triazin-6(5H)-one"),
    ("[H]Sc1nc2cnnnc2[nH]1",         "1H-imidazo[4,5-d][1,2,3]triazin-6(5H)-thione"),
    ("[H]Oc1nnc2nc[nH]c2n1",         "1H-imidazo[5,4-e][1,2,4]triazin-3(1H)-one"),
    ("[H]Sc1nnc2nc[nH]c2n1",         "1H-imidazo[5,4-e][1,2,4]triazin-3(1H)-thione"),
    ("[H]Oc1nc2nncnc2[nH]1",         "1H-imidazo[5,4-e][1,2,4]triazin-6(1H)-one"),
    ("[H]Sc1nc2nncnc2[nH]1",         "1H-imidazo[5,4-e][1,2,4]triazin-6(1H)-thione"),
    ("c1cc2[nH]cnc2nn1",             "5H-imidazo[4,5-c]pyridazine"),
    ("c1nc2cnncc2[nH]1",             "1H-imidazo[4,5-d]pyridazine"),
    ("c1nc2nnncc2[nH]1",             "5H-imidazo[4,5-d][1,2,3]triazine"),
    ("c1nnc2[nH]cnc2n1",             "7H-imidazo[4,5-e][1,2,4]triazine"),
    ("c1nc2cnnnc2[nH]1",             "7H-imidazo[4,5-d][1,2,3]triazine"),
    ("c1nnc2nc[nH]c2n1",             "1H-imidazo[5,4-e][1,2,4]triazine"),
])
def test_phase835(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
