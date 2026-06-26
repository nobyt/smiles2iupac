"""Phase 707: methyl derivatives of oxazolo[x,y-d/e] and thiazolo[x,y-d/e] series
(12 parent compounds).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # oxazolo[4,5-d]pyrimidine (CH at 2,5,7)
    ("c1ncc2ocnc2n1",      "oxazolo[4,5-d]pyrimidine"),
    ("Cc1nc2ncncc2o1",     "2-methyloxazolo[4,5-d]pyrimidine"),
    ("Cc1ncc2ocnc2n1",     "5-methyloxazolo[4,5-d]pyrimidine"),
    ("Cc1ncnc2ncoc12",     "7-methyloxazolo[4,5-d]pyrimidine"),
    # oxazolo[5,4-d]pyrimidine (CH at 2,5,7)
    ("c1ncc2ncoc2n1",      "oxazolo[5,4-d]pyrimidine"),
    ("Cc1nc2cncnc2o1",     "2-methyloxazolo[5,4-d]pyrimidine"),
    ("Cc1ncc2ncoc2n1",     "5-methyloxazolo[5,4-d]pyrimidine"),
    ("Cc1ncnc2ocnc12",     "7-methyloxazolo[5,4-d]pyrimidine"),
    # oxazolo[4,5-d]pyridazine (CH at 2,4,7)
    ("c1nc2cnncc2o1",      "oxazolo[4,5-d]pyridazine"),
    ("Cc1nc2cnncc2o1",     "2-methyloxazolo[4,5-d]pyridazine"),
    ("Cc1nncc2ocnc12",     "4-methyloxazolo[4,5-d]pyridazine"),
    ("Cc1nncc2ncoc12",     "7-methyloxazolo[4,5-d]pyridazine"),
    # oxazolo[4,5-e]pyrazine (CH at 2,5,6)
    ("c1cnc2ocnc2n1",      "oxazolo[4,5-e]pyrazine"),
    ("Cc1nc2nccnc2o1",     "2-methyloxazolo[4,5-e]pyrazine"),
    ("Cc1cnc2ocnc2n1",     "5-methyloxazolo[4,5-e]pyrazine"),
    ("Cc1cnc2ncoc2n1",     "6-methyloxazolo[4,5-e]pyrazine"),
    # oxazolo[4,5-e][1,2,4]triazine (CH at 3,6)
    ("c1nnc2ocnc2n1",      "oxazolo[4,5-e][1,2,4]triazine"),
    ("Cc1nnc2ocnc2n1",     "3-methyloxazolo[4,5-e][1,2,4]triazine"),
    ("Cc1nc2ncnnc2o1",     "6-methyloxazolo[4,5-e][1,2,4]triazine"),
    # oxazolo[5,4-e][1,2,4]triazine (CH at 3,6)
    ("c1nnc2ncoc2n1",      "oxazolo[5,4-e][1,2,4]triazine"),
    ("Cc1nnc2ncoc2n1",     "3-methyloxazolo[5,4-e][1,2,4]triazine"),
    ("Cc1nc2nncnc2o1",     "6-methyloxazolo[5,4-e][1,2,4]triazine"),
    # thiazolo[4,5-d]pyrimidine (CH at 2,5,7)
    ("c1ncc2scnc2n1",      "thiazolo[4,5-d]pyrimidine"),
    ("Cc1nc2ncncc2s1",     "2-methylthiazolo[4,5-d]pyrimidine"),
    ("Cc1ncc2scnc2n1",     "5-methylthiazolo[4,5-d]pyrimidine"),
    ("Cc1ncnc2ncsc12",     "7-methylthiazolo[4,5-d]pyrimidine"),
    # thiazolo[5,4-d]pyrimidine (CH at 2,5,7)
    ("c1ncc2ncsc2n1",      "thiazolo[5,4-d]pyrimidine"),
    ("Cc1nc2cncnc2s1",     "2-methylthiazolo[5,4-d]pyrimidine"),
    ("Cc1ncc2ncsc2n1",     "5-methylthiazolo[5,4-d]pyrimidine"),
    ("Cc1ncnc2scnc12",     "7-methylthiazolo[5,4-d]pyrimidine"),
    # thiazolo[4,5-d]pyridazine (CH at 2,4,7)
    ("c1nc2cnncc2s1",      "thiazolo[4,5-d]pyridazine"),
    ("Cc1nc2cnncc2s1",     "2-methylthiazolo[4,5-d]pyridazine"),
    ("Cc1nncc2scnc12",     "4-methylthiazolo[4,5-d]pyridazine"),
    ("Cc1nncc2ncsc12",     "7-methylthiazolo[4,5-d]pyridazine"),
    # thiazolo[4,5-e]pyrazine (CH at 2,5,6)
    ("c1cnc2scnc2n1",      "thiazolo[4,5-e]pyrazine"),
    ("Cc1nc2nccnc2s1",     "2-methylthiazolo[4,5-e]pyrazine"),
    ("Cc1cnc2scnc2n1",     "5-methylthiazolo[4,5-e]pyrazine"),
    ("Cc1cnc2ncsc2n1",     "6-methylthiazolo[4,5-e]pyrazine"),
    # thiazolo[4,5-e][1,2,4]triazine (CH at 3,6)
    ("c1nnc2scnc2n1",      "thiazolo[4,5-e][1,2,4]triazine"),
    ("Cc1nnc2scnc2n1",     "3-methylthiazolo[4,5-e][1,2,4]triazine"),
    ("Cc1nc2ncnnc2s1",     "6-methylthiazolo[4,5-e][1,2,4]triazine"),
    # thiazolo[5,4-e][1,2,4]triazine (CH at 3,6)
    ("c1nnc2ncsc2n1",      "thiazolo[5,4-e][1,2,4]triazine"),
    ("Cc1nnc2ncsc2n1",     "3-methylthiazolo[5,4-e][1,2,4]triazine"),
    ("Cc1nc2nncnc2s1",     "6-methylthiazolo[5,4-e][1,2,4]triazine"),
])
def test_phase707(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
