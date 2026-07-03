"""Phase 834: 1H-pyrazolo[5,4-c/d] α-ol/thiol → tautomers."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("[H]Oc1n[nH]c2nnccc12",         "1H-pyrazolo[5,4-c]pyridazin-3(1H)-one"),
    ("[H]Sc1n[nH]c2nnccc12",         "1H-pyrazolo[5,4-c]pyridazin-3(1H)-thione"),
    ("[H]Oc1cc2cn[nH]c2nn1",         "1H-pyrazolo[5,4-c]pyridazin-5(3H)-one"),
    ("[H]Sc1cc2cn[nH]c2nn1",         "1H-pyrazolo[5,4-c]pyridazin-5(3H)-thione"),
    ("[H]Oc1nnnc2[nH]ncc12",         "1H-pyrazolo[5,4-d][1,2,3]triazin-4(4H)-one"),
    ("[H]Sc1nnnc2[nH]ncc12",         "1H-pyrazolo[5,4-d][1,2,3]triazin-4(4H)-thione"),
    ("[H]Oc1n[nH]c2nnncc12",         "1H-pyrazolo[5,4-d][1,2,3]triazin-5(1H)-one"),
    ("[H]Sc1n[nH]c2nnncc12",         "1H-pyrazolo[5,4-d][1,2,3]triazin-5(1H)-thione"),
    ("c1cc2cn[nH]c2nn1",             "1H-pyrazolo[5,4-c]pyridazine"),
    ("c1nnnc2[nH]ncc12",             "1H-pyrazolo[5,4-d][1,2,3]triazine"),
])
def test_phase834(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
