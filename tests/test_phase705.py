"""Phase 705: methyl derivatives of furo[x,y-d/e] bicyclic series
(7 parent compounds: pyrimidine, pyridazine, pyrazine, triazine partners).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # furo[2,3-d]pyrimidine (CH at 2,4,5,6)
    ("c1ncc2ccoc2n1",      "furo[2,3-d]pyrimidine"),
    ("Cc1ncc2ccoc2n1",     "2-methylfuro[2,3-d]pyrimidine"),
    ("Cc1ncnc2occc12",     "4-methylfuro[2,3-d]pyrimidine"),
    ("Cc1coc2ncncc12",     "5-methylfuro[2,3-d]pyrimidine"),
    ("Cc1cc2cncnc2o1",     "6-methylfuro[2,3-d]pyrimidine"),
    # furo[3,2-d]pyrimidine (CH at 2,4,6,7)
    ("c1ncc2occc2n1",      "furo[3,2-d]pyrimidine"),
    ("Cc1ncc2occc2n1",     "2-methylfuro[3,2-d]pyrimidine"),
    ("Cc1ncnc2ccoc12",     "4-methylfuro[3,2-d]pyrimidine"),
    ("Cc1cc2ncncc2o1",     "6-methylfuro[3,2-d]pyrimidine"),
    ("Cc1coc2cncnc12",     "7-methylfuro[3,2-d]pyrimidine"),
    # furo[3,2-d]pyridazine (CH at 2,3,4,7)
    ("c1cc2cnncc2o1",      "furo[3,2-d]pyridazine"),
    ("Cc1cc2cnncc2o1",     "2-methylfuro[3,2-d]pyridazine"),
    ("Cc1coc2cnncc12",     "3-methylfuro[3,2-d]pyridazine"),
    ("Cc1nncc2occc12",     "4-methylfuro[3,2-d]pyridazine"),
    ("Cc1nncc2ccoc12",     "7-methylfuro[3,2-d]pyridazine"),
    # furo[2,3-e]pyrazine (CH at 2,3,6,7)
    ("c1cnc2occc2n1",      "furo[2,3-e]pyrazine"),
    ("Cc1cnc2occc2n1",     "2-methylfuro[2,3-e]pyrazine"),
    ("Cc1cnc2ccoc2n1",     "3-methylfuro[2,3-e]pyrazine"),
    ("Cc1cc2nccnc2o1",     "6-methylfuro[2,3-e]pyrazine"),
    ("Cc1coc2nccnc12",     "7-methylfuro[2,3-e]pyrazine"),
    # furo[3,4-d]pyrimidine (CH at 2,4,5,7)
    ("c1ncc2cocc2n1",      "furo[3,4-d]pyrimidine"),
    ("Cc1ncc2cocc2n1",     "2-methylfuro[3,4-d]pyrimidine"),
    ("Cc1ncnc2cocc12",     "4-methylfuro[3,4-d]pyrimidine"),
    ("Cc1occ2ncncc12",     "5-methylfuro[3,4-d]pyrimidine"),
    ("Cc1occ2cncnc12",     "7-methylfuro[3,4-d]pyrimidine"),
    # furo[2,3-d][1,2,3]triazine (CH at 4,5,6)
    ("c1cc2cnnnc2o1",      "furo[2,3-d][1,2,3]triazine"),
    ("Cc1nnnc2occc12",     "4-methylfuro[2,3-d][1,2,3]triazine"),
    ("Cc1coc2nnncc12",     "5-methylfuro[2,3-d][1,2,3]triazine"),
    ("Cc1cc2cnnnc2o1",     "6-methylfuro[2,3-d][1,2,3]triazine"),
    # furo[2,3-e][1,2,4]triazine (CH at 3,6,7)
    ("c1nnc2ccoc2n1",      "furo[2,3-e][1,2,4]triazine"),
    ("Cc1nnc2ccoc2n1",     "3-methylfuro[2,3-e][1,2,4]triazine"),
    ("Cc1cc2nncnc2o1",     "6-methylfuro[2,3-e][1,2,4]triazine"),
    ("Cc1coc2ncnnc12",     "7-methylfuro[2,3-e][1,2,4]triazine"),
])
def test_phase705(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
