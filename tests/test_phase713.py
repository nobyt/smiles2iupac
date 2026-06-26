"""Phase 713: methyl derivatives of furo[x,y-b/c] series and remaining furo d/e
(15 parent compounds: b/c-fused with pyridine/pyridazine, and
remaining 3,2- and 3,4-series d/e isomers).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # furo[2,3-b]pyridine (CH at 2,3,4,5,6)
    ("c1cnc2occc2c1",       "furo[2,3-b]pyridine"),
    ("Cc1cc2cccnc2o1",      "2-methylfuro[2,3-b]pyridine"),
    ("Cc1coc2ncccc12",      "3-methylfuro[2,3-b]pyridine"),
    ("Cc1ccnc2occc12",      "4-methylfuro[2,3-b]pyridine"),
    ("Cc1cnc2occc2c1",      "5-methylfuro[2,3-b]pyridine"),
    ("Cc1ccc2ccoc2n1",      "6-methylfuro[2,3-b]pyridine"),
    # furo[2,3-c]pyridine (CH at 2,3,4,5,7)
    ("c1cc2ccoc2cn1",       "furo[2,3-c]pyridine"),
    ("Cc1cc2ccncc2o1",      "2-methylfuro[2,3-c]pyridine"),
    ("Cc1coc2cnccc12",      "3-methylfuro[2,3-c]pyridine"),
    ("Cc1cncc2occc12",      "4-methylfuro[2,3-c]pyridine"),
    ("Cc1cc2ccoc2cn1",      "5-methylfuro[2,3-c]pyridine"),
    ("Cc1nccc2ccoc12",      "7-methylfuro[2,3-c]pyridine"),
    # furo[2,3-c]pyridazine (CH at 3,4,5,6)
    ("c1cc2ccoc2nn1",       "furo[2,3-c]pyridazine"),
    ("Cc1cc2ccoc2nn1",      "3-methylfuro[2,3-c]pyridazine"),
    ("Cc1cnnc2occc12",      "4-methylfuro[2,3-c]pyridazine"),
    ("Cc1coc2nnccc12",      "5-methylfuro[2,3-c]pyridazine"),
    ("Cc1cc2ccnnc2o1",      "6-methylfuro[2,3-c]pyridazine"),
    # furo[2,3-e]pyridazine (CH at 3,4,6,7)
    ("c1cc2occc2nn1",       "furo[2,3-e]pyridazine"),
    ("Cc1cc2occc2nn1",      "3-methylfuro[2,3-e]pyridazine"),
    ("Cc1cnnc2ccoc12",      "4-methylfuro[2,3-e]pyridazine"),
    ("Cc1cc2nnccc2o1",      "6-methylfuro[2,3-e]pyridazine"),
    ("Cc1coc2ccnnc12",      "7-methylfuro[2,3-e]pyridazine"),
    # furo[3,2-b]pyridine (CH at 2,3,5,6,7)
    ("c1cnc2ccoc2c1",       "furo[3,2-b]pyridine"),
    ("Cc1cc2ncccc2o1",      "2-methylfuro[3,2-b]pyridine"),
    ("Cc1coc2cccnc12",      "3-methylfuro[3,2-b]pyridine"),
    ("Cc1ccc2occc2n1",      "5-methylfuro[3,2-b]pyridine"),
    ("Cc1cnc2ccoc2c1",      "6-methylfuro[3,2-b]pyridine"),
    ("Cc1ccnc2ccoc12",      "7-methylfuro[3,2-b]pyridine"),
    # furo[3,2-c]pyridine (CH at 2,3,4,6,7)
    ("c1cc2occc2cn1",       "furo[3,2-c]pyridine"),
    ("Cc1cc2cnccc2o1",      "2-methylfuro[3,2-c]pyridine"),
    ("Cc1coc2ccncc12",      "3-methylfuro[3,2-c]pyridine"),
    ("Cc1nccc2occc12",      "4-methylfuro[3,2-c]pyridine"),
    ("Cc1cc2occc2cn1",      "6-methylfuro[3,2-c]pyridine"),
    ("Cc1cncc2ccoc12",      "7-methylfuro[3,2-c]pyridine"),
    # furo[3,2-d][1,2,3]triazine (CH at 4,6,7)
    ("c1cc2nnncc2o1",       "furo[3,2-d][1,2,3]triazine"),
    ("Cc1nnnc2ccoc12",      "4-methylfuro[3,2-d][1,2,3]triazine"),
    ("Cc1cc2nnncc2o1",      "6-methylfuro[3,2-d][1,2,3]triazine"),
    ("Cc1coc2cnnnc12",      "7-methylfuro[3,2-d][1,2,3]triazine"),
    # furo[3,2-e][1,2,4]triazine (CH at 3,5,6)
    ("c1nnc2occc2n1",       "furo[3,2-e][1,2,4]triazine"),
    ("Cc1nnc2occc2n1",      "3-methylfuro[3,2-e][1,2,4]triazine"),
    ("Cc1coc2nncnc12",      "5-methylfuro[3,2-e][1,2,4]triazine"),
    ("Cc1cc2ncnnc2o1",      "6-methylfuro[3,2-e][1,2,4]triazine"),
    # furo[3,4-b]pyridine (CH at 2,3,4,5,7)
    ("c1cnc2cocc2c1",       "furo[3,4-b]pyridine"),
    ("Cc1ccc2cocc2n1",      "2-methylfuro[3,4-b]pyridine"),
    ("Cc1cnc2cocc2c1",      "3-methylfuro[3,4-b]pyridine"),
    ("Cc1ccnc2cocc12",      "4-methylfuro[3,4-b]pyridine"),
    ("Cc1occ2ncccc12",      "5-methylfuro[3,4-b]pyridine"),
    ("Cc1occ2cccnc12",      "7-methylfuro[3,4-b]pyridine"),
    # furo[3,4-c]pyridine (CH at 1,3,4,6,7)
    ("c1cc2cocc2cn1",       "furo[3,4-c]pyridine"),
    ("Cc1occ2cnccc12",      "1-methylfuro[3,4-c]pyridine"),
    ("Cc1occ2ccncc12",      "3-methylfuro[3,4-c]pyridine"),
    ("Cc1nccc2cocc12",      "4-methylfuro[3,4-c]pyridine"),
    ("Cc1cc2cocc2cn1",      "6-methylfuro[3,4-c]pyridine"),
    ("Cc1cncc2cocc12",      "7-methylfuro[3,4-c]pyridine"),
    # furo[3,4-c]pyridazine (CH at 3,4,5,7)
    ("c1cc2cocc2nn1",       "furo[3,4-c]pyridazine"),
    ("Cc1cc2cocc2nn1",      "3-methylfuro[3,4-c]pyridazine"),
    ("Cc1cnnc2cocc12",      "4-methylfuro[3,4-c]pyridazine"),
    ("Cc1occ2nnccc12",      "5-methylfuro[3,4-c]pyridazine"),
    ("Cc1occ2ccnnc12",      "7-methylfuro[3,4-c]pyridazine"),
    # furo[3,4-d][1,2,3]triazine (CH at 4,5,7)
    ("c1nnnc2cocc12",       "furo[3,4-d][1,2,3]triazine"),
    ("Cc1nnnc2cocc12",      "4-methylfuro[3,4-d][1,2,3]triazine"),
    ("Cc1occ2nnncc12",      "5-methylfuro[3,4-d][1,2,3]triazine"),
    ("Cc1occ2cnnnc12",      "7-methylfuro[3,4-d][1,2,3]triazine"),
    # furo[3,4-d]pyridazine (CH at 1,5)
    ("c1nncc2cocc12",       "furo[3,4-d]pyridazine"),
    ("Cc1nncc2cocc12",      "1-methylfuro[3,4-d]pyridazine"),
    ("Cc1occ2cnncc12",      "5-methylfuro[3,4-d]pyridazine"),
    # furo[3,4-e][1,2,4]triazine (CH at 3,5,7)
    ("c1nnc2cocc2n1",       "furo[3,4-e][1,2,4]triazine"),
    ("Cc1nnc2cocc2n1",      "3-methylfuro[3,4-e][1,2,4]triazine"),
    ("Cc1occ2nncnc12",      "5-methylfuro[3,4-e][1,2,4]triazine"),
    ("Cc1occ2ncnnc12",      "7-methylfuro[3,4-e][1,2,4]triazine"),
    # furo[3,4-e]pyrazine (CH at 2,5)
    ("c1cnc2cocc2n1",       "furo[3,4-e]pyrazine"),
    ("Cc1cnc2cocc2n1",      "2-methylfuro[3,4-e]pyrazine"),
    ("Cc1occ2nccnc12",      "5-methylfuro[3,4-e]pyrazine"),
])
def test_phase713(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
