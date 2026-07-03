"""Phase 830: 1H-pyrrolo[3,2-c/d/e] α-ol/thiol → tautomers."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("[H]Oc1cc2nccc-2[nH]n1",    "1H-pyrrolo[3,2-c]pyridazin-3(1H)-one"),
    ("[H]Sc1cc2nccc-2[nH]n1",    "1H-pyrrolo[3,2-c]pyridazin-3(1H)-thione"),
    ("[H]Oc1cc2[nH]nccc-2n1",    "1H-pyrrolo[3,2-c]pyridazin-6(1H)-one"),
    ("[H]Sc1cc2[nH]nccc-2n1",    "1H-pyrrolo[3,2-c]pyridazin-6(1H)-thione"),
    ("[H]Oc1ncc2nccc-2[nH]1",    "1H-pyrrolo[3,2-d]pyrimidin-2(1H)-one"),
    ("[H]Sc1ncc2nccc-2[nH]1",    "1H-pyrrolo[3,2-d]pyrimidin-2(1H)-thione"),
    ("[H]Oc1nc[nH]c2ccnc1-2",    "1H-pyrrolo[3,2-d]pyrimidin-4(4H)-one"),
    ("[H]Sc1nc[nH]c2ccnc1-2",    "1H-pyrrolo[3,2-d]pyrimidin-4(4H)-thione"),
    ("[H]Oc1cc2[nH]cncc-2n1",    "1H-pyrrolo[3,2-d]pyrimidin-6(1H)-one"),
    ("[H]Sc1cc2[nH]cncc-2n1",    "1H-pyrrolo[3,2-d]pyrimidin-6(1H)-thione"),
    ("[H]Oc1nn[nH]c2ccnc1-2",    "1H-pyrrolo[3,2-d][1,2,3]triazin-4(4H)-one"),
    ("[H]Sc1nn[nH]c2ccnc1-2",    "1H-pyrrolo[3,2-d][1,2,3]triazin-4(4H)-thione"),
    ("[H]Oc1cc2[nH]nncc-2n1",    "1H-pyrrolo[3,2-d][1,2,3]triazin-6(1H)-one"),
    ("[H]Sc1cc2[nH]nncc-2n1",    "1H-pyrrolo[3,2-d][1,2,3]triazin-6(1H)-thione"),
    ("[H]Oc1nnc2[nH]ccc2n1",      "1H-pyrrolo[3,2-e][1,2,4]triazin-3(1H)-one"),
    ("[H]Sc1nnc2[nH]ccc2n1",      "1H-pyrrolo[3,2-e][1,2,4]triazin-3(1H)-thione"),
    ("[H]Oc1cc2ncnnc2[nH]1",      "1H-pyrrolo[3,2-e][1,2,4]triazin-6(1H)-one"),
    ("[H]Sc1cc2ncnnc2[nH]1",      "1H-pyrrolo[3,2-e][1,2,4]triazin-6(1H)-thione"),
    ("c1cc2nccc-2[nH]n1",          "1H-pyrrolo[3,2-c]pyridazine"),
    ("c1cc2[nH]cncc-2n1",          "1H-pyrrolo[3,2-d]pyrimidine"),
    ("c1cc2[nH]nncc-2n1",          "1H-pyrrolo[3,2-d][1,2,3]triazine"),
    ("c1nnc2[nH]ccc2n1",           "1H-pyrrolo[3,2-e][1,2,4]triazine"),
])
def test_phase830(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
