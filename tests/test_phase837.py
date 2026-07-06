"""Phase 837: 2H-[1,2,3]triazolo[4,5-c/d/e] α-ol/thiol → tautomers."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("[H]Oc1cc2n[nH]nc2nn1",         "2H-[1,2,3]triazolo[4,5-c]pyridazin-6(5H)-one"),
    ("[H]Sc1cc2n[nH]nc2nn1",         "2H-[1,2,3]triazolo[4,5-c]pyridazin-6(5H)-thione"),
    ("[H]Oc1nccc2n[nH]nc12",         "2H-[1,2,3]triazolo[4,5-c]pyridin-4(3H)-one"),
    ("[H]Sc1nccc2n[nH]nc12",         "2H-[1,2,3]triazolo[4,5-c]pyridin-4(3H)-thione"),
    ("[H]Oc1cc2n[nH]nc2cn1",         "2H-[1,2,3]triazolo[4,5-c]pyridin-6(5H)-one"),
    ("[H]Sc1cc2n[nH]nc2cn1",         "2H-[1,2,3]triazolo[4,5-c]pyridin-6(5H)-thione"),
    ("[H]Oc1nnnc2n[nH]nc12",         "2H-[1,2,3]triazolo[4,5-d][1,2,3]triazin-7(3H)-one"),
    ("[H]Sc1nnnc2n[nH]nc12",         "2H-[1,2,3]triazolo[4,5-d][1,2,3]triazin-7(3H)-thione"),
    ("[H]Oc1ncc2n[nH]nc2n1",         "2H-[1,2,3]triazolo[4,5-d]pyrimidin-5(3H)-one"),
    ("[H]Sc1ncc2n[nH]nc2n1",         "2H-[1,2,3]triazolo[4,5-d]pyrimidin-5(3H)-thione"),
    ("[H]Oc1ncnc2n[nH]nc12",         "2H-[1,2,3]triazolo[4,5-d]pyrimidin-7(3H)-one"),
    ("[H]Sc1ncnc2n[nH]nc12",         "2H-[1,2,3]triazolo[4,5-d]pyrimidin-7(3H)-thione"),
    ("[H]Oc1nnc2n[nH]nc2n1",         "2H-[1,2,3]triazolo[4,5-e][1,2,4]triazin-6(5H)-one"),
    ("[H]Sc1nnc2n[nH]nc2n1",         "2H-[1,2,3]triazolo[4,5-e][1,2,4]triazin-6(5H)-thione"),
    ("[H]Oc1cnc2n[nH]nc2n1",         "2H-[1,2,3]triazolo[4,5-e]pyrazin-5(3H)-one"),
    ("[H]Sc1cnc2n[nH]nc2n1",         "2H-[1,2,3]triazolo[4,5-e]pyrazin-5(3H)-thione"),
    ("c1cc2n[nH]nc2nn1",             "2H-[1,2,3]triazolo[4,5-c]pyridazine"),
    ("c1cc2n[nH]nc2cn1",             "2H-[1,2,3]triazolo[4,5-c]pyridine"),
    ("c1nnnc2n[nH]nc12",             "2H-[1,2,3]triazolo[4,5-d][1,2,3]triazine"),
    ("c1ncc2n[nH]nc2n1",             "2H-[1,2,3]triazolo[4,5-d]pyrimidine"),
    ("c1nnc2n[nH]nc2n1",             "2H-[1,2,3]triazolo[4,5-e][1,2,4]triazine"),
    ("c1cnc2n[nH]nc2n1",             "2H-[1,2,3]triazolo[4,5-e]pyrazine"),
])
def test_phase837(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
