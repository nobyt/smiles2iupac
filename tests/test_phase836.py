"""Phase 836: 1H-[1,2,3]triazolo[4,5-c/d/e] + 1H-[1,2,3]triazolo[5,4-c/d/e] α-ol/thiol → tautomers."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("[H]Oc1cc2[nH]nnc2nn1",         "1H-[1,2,3]triazolo[4,5-c]pyridazin-6(1H)-one"),
    ("[H]Sc1cc2[nH]nnc2nn1",         "1H-[1,2,3]triazolo[4,5-c]pyridazin-6(1H)-thione"),
    ("[H]Oc1nnnc2nn[nH]c12",         "1H-[1,2,3]triazolo[4,5-d][1,2,3]triazin-7(3H)-one"),
    ("[H]Sc1nnnc2nn[nH]c12",         "1H-[1,2,3]triazolo[4,5-d][1,2,3]triazin-7(3H)-thione"),
    ("[H]Oc1ncc2[nH]nnc2n1",         "1H-[1,2,3]triazolo[4,5-d]pyrimidin-5(3H)-one"),
    ("[H]Sc1ncc2[nH]nnc2n1",         "1H-[1,2,3]triazolo[4,5-d]pyrimidin-5(3H)-thione"),
    ("[H]Oc1ncnc2nn[nH]c12",         "1H-[1,2,3]triazolo[4,5-d]pyrimidin-7(3H)-one"),
    ("[H]Sc1ncnc2nn[nH]c12",         "1H-[1,2,3]triazolo[4,5-d]pyrimidin-7(3H)-thione"),
    ("[H]Oc1nnc2[nH]nnc2n1",         "1H-[1,2,3]triazolo[4,5-e][1,2,4]triazin-6(1H)-one"),
    ("[H]Sc1nnc2[nH]nnc2n1",         "1H-[1,2,3]triazolo[4,5-e][1,2,4]triazin-6(1H)-thione"),
    ("[H]Oc1cc2nn[nH]c2nn1",         "1H-[1,2,3]triazolo[5,4-c]pyridazin-6(1H)-one"),
    ("[H]Sc1cc2nn[nH]c2nn1",         "1H-[1,2,3]triazolo[5,4-c]pyridazin-6(1H)-thione"),
    ("[H]Oc1nccc2nn[nH]c12",         "1H-[1,2,3]triazolo[5,4-c]pyridin-4(1H)-one"),
    ("[H]Sc1nccc2nn[nH]c12",         "1H-[1,2,3]triazolo[5,4-c]pyridin-4(1H)-thione"),
    ("[H]Oc1cc2nn[nH]c2cn1",         "1H-[1,2,3]triazolo[5,4-c]pyridin-6(1H)-one"),
    ("[H]Sc1cc2nn[nH]c2cn1",         "1H-[1,2,3]triazolo[5,4-c]pyridin-6(1H)-thione"),
    ("[H]Oc1nnnc2[nH]nnc12",         "1H-[1,2,3]triazolo[5,4-d][1,2,3]triazin-7(3H)-one"),
    ("[H]Sc1nnnc2[nH]nnc12",         "1H-[1,2,3]triazolo[5,4-d][1,2,3]triazin-7(3H)-thione"),
    ("[H]Oc1ncc2nn[nH]c2n1",         "1H-[1,2,3]triazolo[5,4-d]pyrimidin-5(3H)-one"),
    ("[H]Sc1ncc2nn[nH]c2n1",         "1H-[1,2,3]triazolo[5,4-d]pyrimidin-5(3H)-thione"),
    ("[H]Oc1ncnc2[nH]nnc12",         "1H-[1,2,3]triazolo[5,4-d]pyrimidin-7(3H)-one"),
    ("[H]Sc1ncnc2[nH]nnc12",         "1H-[1,2,3]triazolo[5,4-d]pyrimidin-7(3H)-thione"),
    ("[H]Oc1nnc2nn[nH]c2n1",         "1H-[1,2,3]triazolo[5,4-e][1,2,4]triazin-6(1H)-one"),
    ("[H]Sc1nnc2nn[nH]c2n1",         "1H-[1,2,3]triazolo[5,4-e][1,2,4]triazin-6(1H)-thione"),
    ("c1cc2[nH]nnc2nn1",             "1H-[1,2,3]triazolo[4,5-c]pyridazine"),
    ("c1nnnc2nn[nH]c12",             "1H-[1,2,3]triazolo[4,5-d][1,2,3]triazine"),
    ("c1ncc2[nH]nnc2n1",             "1H-[1,2,3]triazolo[4,5-d]pyrimidine"),
    ("c1nnc2[nH]nnc2n1",             "1H-[1,2,3]triazolo[4,5-e][1,2,4]triazine"),
    ("c1cc2nn[nH]c2nn1",             "1H-[1,2,3]triazolo[5,4-c]pyridazine"),
    ("c1cc2nn[nH]c2cn1",             "1H-[1,2,3]triazolo[5,4-c]pyridine"),
    ("c1nnnc2[nH]nnc12",             "1H-[1,2,3]triazolo[5,4-d][1,2,3]triazine"),
    ("c1ncc2nn[nH]c2n1",             "1H-[1,2,3]triazolo[5,4-d]pyrimidine"),
    ("c1nnc2nn[nH]c2n1",             "1H-[1,2,3]triazolo[5,4-e][1,2,4]triazine"),
])
def test_phase836(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
