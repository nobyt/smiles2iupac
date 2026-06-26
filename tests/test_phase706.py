"""Phase 706: methyl derivatives of thieno[x,y-d/e] bicyclic series
(10 parent compounds: pyrimidine, pyridazine, pyrazine, triazine partners).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[2,3-d]pyrimidine (CH at 2,4,5,6)
    ("c1ncc2ccsc2n1",      "thieno[2,3-d]pyrimidine"),
    ("Cc1ncc2ccsc2n1",     "2-methylthieno[2,3-d]pyrimidine"),
    ("Cc1ncnc2sccc12",     "4-methylthieno[2,3-d]pyrimidine"),
    ("Cc1csc2ncncc12",     "5-methylthieno[2,3-d]pyrimidine"),
    ("Cc1cc2cncnc2s1",     "6-methylthieno[2,3-d]pyrimidine"),
    # thieno[3,2-d]pyrimidine (CH at 2,4,6,7)
    ("c1ncc2sccc2n1",      "thieno[3,2-d]pyrimidine"),
    ("Cc1ncc2sccc2n1",     "2-methylthieno[3,2-d]pyrimidine"),
    ("Cc1ncnc2ccsc12",     "4-methylthieno[3,2-d]pyrimidine"),
    ("Cc1cc2ncncc2s1",     "6-methylthieno[3,2-d]pyrimidine"),
    ("Cc1csc2cncnc12",     "7-methylthieno[3,2-d]pyrimidine"),
    # thieno[3,2-d]pyridazine (CH at 2,3,4,7)
    ("c1cc2cnncc2s1",      "thieno[3,2-d]pyridazine"),
    ("Cc1cc2cnncc2s1",     "2-methylthieno[3,2-d]pyridazine"),
    ("Cc1csc2cnncc12",     "3-methylthieno[3,2-d]pyridazine"),
    ("Cc1nncc2sccc12",     "4-methylthieno[3,2-d]pyridazine"),
    ("Cc1nncc2ccsc12",     "7-methylthieno[3,2-d]pyridazine"),
    # thieno[2,3-e]pyrazine (CH at 2,3,6,7)
    ("c1cnc2sccc2n1",      "thieno[2,3-e]pyrazine"),
    ("Cc1cnc2sccc2n1",     "2-methylthieno[2,3-e]pyrazine"),
    ("Cc1cnc2ccsc2n1",     "3-methylthieno[2,3-e]pyrazine"),
    ("Cc1cc2nccnc2s1",     "6-methylthieno[2,3-e]pyrazine"),
    ("Cc1csc2nccnc12",     "7-methylthieno[2,3-e]pyrazine"),
    # thieno[3,4-d]pyrimidine (CH at 2,4,5,7)
    ("c1ncc2cscc2n1",      "thieno[3,4-d]pyrimidine"),
    ("Cc1ncc2cscc2n1",     "2-methylthieno[3,4-d]pyrimidine"),
    ("Cc1ncnc2cscc12",     "4-methylthieno[3,4-d]pyrimidine"),
    ("Cc1scc2ncncc12",     "5-methylthieno[3,4-d]pyrimidine"),
    ("Cc1scc2cncnc12",     "7-methylthieno[3,4-d]pyrimidine"),
    # thieno[2,3-d][1,2,3]triazine (CH at 4,5,6)
    ("c1cc2cnnnc2s1",      "thieno[2,3-d][1,2,3]triazine"),
    ("Cc1nnnc2sccc12",     "4-methylthieno[2,3-d][1,2,3]triazine"),
    ("Cc1csc2nnncc12",     "5-methylthieno[2,3-d][1,2,3]triazine"),
    ("Cc1cc2cnnnc2s1",     "6-methylthieno[2,3-d][1,2,3]triazine"),
    # thieno[2,3-e][1,2,4]triazine (CH at 3,6,7)
    ("c1nnc2ccsc2n1",      "thieno[2,3-e][1,2,4]triazine"),
    ("Cc1nnc2ccsc2n1",     "3-methylthieno[2,3-e][1,2,4]triazine"),
    ("Cc1cc2nncnc2s1",     "6-methylthieno[2,3-e][1,2,4]triazine"),
    ("Cc1csc2ncnnc12",     "7-methylthieno[2,3-e][1,2,4]triazine"),
    # thieno[3,2-e][1,2,4]triazine (CH at 3,5,6)
    ("c1nnc2sccc2n1",      "thieno[3,2-e][1,2,4]triazine"),
    ("Cc1nnc2sccc2n1",     "3-methylthieno[3,2-e][1,2,4]triazine"),
    ("Cc1csc2nncnc12",     "5-methylthieno[3,2-e][1,2,4]triazine"),
    ("Cc1cc2ncnnc2s1",     "6-methylthieno[3,2-e][1,2,4]triazine"),
    # thieno[3,4-e][1,2,4]triazine (CH at 3,5,7)
    ("c1nnc2cscc2n1",      "thieno[3,4-e][1,2,4]triazine"),
    ("Cc1nnc2cscc2n1",     "3-methylthieno[3,4-e][1,2,4]triazine"),
    ("Cc1scc2nncnc12",     "5-methylthieno[3,4-e][1,2,4]triazine"),
    ("Cc1scc2ncnnc12",     "7-methylthieno[3,4-e][1,2,4]triazine"),
    # thieno[3,4-d][1,2,3]triazine (CH at 4,5,7)
    ("c1nnnc2cscc12",      "thieno[3,4-d][1,2,3]triazine"),
    ("Cc1nnnc2cscc12",     "4-methylthieno[3,4-d][1,2,3]triazine"),
    ("Cc1scc2nnncc12",     "5-methylthieno[3,4-d][1,2,3]triazine"),
    ("Cc1scc2cnnnc12",     "7-methylthieno[3,4-d][1,2,3]triazine"),
])
def test_phase706(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
