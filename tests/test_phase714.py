"""Phase 714: methyl derivatives of thieno[x,y-b/c] series and remaining thieno d/e
(17 parent compounds: b/c-fused with pyridine/pyridazine,
remaining 3,2- and 3,4-series d/e isomers, and thienothiophenes).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[2,3-b]pyridine (CH at 2,3,4,5,6)
    ("c1cnc2sccc2c1",       "thieno[2,3-b]pyridine"),
    ("Cc1cc2cccnc2s1",      "2-methylthieno[2,3-b]pyridine"),
    ("Cc1csc2ncccc12",      "3-methylthieno[2,3-b]pyridine"),
    ("Cc1ccnc2sccc12",      "4-methylthieno[2,3-b]pyridine"),
    ("Cc1cnc2sccc2c1",      "5-methylthieno[2,3-b]pyridine"),
    ("Cc1ccc2ccsc2n1",      "6-methylthieno[2,3-b]pyridine"),
    # thieno[2,3-c]pyridine (CH at 2,3,4,5,7)
    ("c1cc2ccsc2cn1",       "thieno[2,3-c]pyridine"),
    ("Cc1cc2ccncc2s1",      "2-methylthieno[2,3-c]pyridine"),
    ("Cc1csc2cnccc12",      "3-methylthieno[2,3-c]pyridine"),
    ("Cc1cncc2sccc12",      "4-methylthieno[2,3-c]pyridine"),
    ("Cc1cc2ccsc2cn1",      "5-methylthieno[2,3-c]pyridine"),
    ("Cc1nccc2ccsc12",      "7-methylthieno[2,3-c]pyridine"),
    # thieno[2,3-c]pyridazine (CH at 3,4,5,6)
    ("c1cc2ccsc2nn1",       "thieno[2,3-c]pyridazine"),
    ("Cc1cc2ccsc2nn1",      "3-methylthieno[2,3-c]pyridazine"),
    ("Cc1cnnc2sccc12",      "4-methylthieno[2,3-c]pyridazine"),
    ("Cc1csc2nnccc12",      "5-methylthieno[2,3-c]pyridazine"),
    ("Cc1cc2ccnnc2s1",      "6-methylthieno[2,3-c]pyridazine"),
    # thieno[3,2-b]pyridine (CH at 2,3,5,6,7)
    ("c1cnc2ccsc2c1",       "thieno[3,2-b]pyridine"),
    ("Cc1cc2ncccc2s1",      "2-methylthieno[3,2-b]pyridine"),
    ("Cc1csc2cccnc12",      "3-methylthieno[3,2-b]pyridine"),
    ("Cc1ccc2sccc2n1",      "5-methylthieno[3,2-b]pyridine"),
    ("Cc1cnc2ccsc2c1",      "6-methylthieno[3,2-b]pyridine"),
    ("Cc1ccnc2ccsc12",      "7-methylthieno[3,2-b]pyridine"),
    # thieno[3,2-c]pyridine (CH at 2,3,4,6,7)
    ("c1cc2sccc2cn1",       "thieno[3,2-c]pyridine"),
    ("Cc1cc2cnccc2s1",      "2-methylthieno[3,2-c]pyridine"),
    ("Cc1csc2ccncc12",      "3-methylthieno[3,2-c]pyridine"),
    ("Cc1nccc2sccc12",      "4-methylthieno[3,2-c]pyridine"),
    ("Cc1cc2sccc2cn1",      "6-methylthieno[3,2-c]pyridine"),
    ("Cc1cncc2ccsc12",      "7-methylthieno[3,2-c]pyridine"),
    # thieno[3,2-c]pyridazine (CH at 3,4,6,7)
    ("c1cc2sccc2nn1",       "thieno[3,2-c]pyridazine"),
    ("Cc1cc2sccc2nn1",      "3-methylthieno[3,2-c]pyridazine"),
    ("Cc1cnnc2ccsc12",      "4-methylthieno[3,2-c]pyridazine"),
    ("Cc1cc2nnccc2s1",      "6-methylthieno[3,2-c]pyridazine"),
    ("Cc1csc2ccnnc12",      "7-methylthieno[3,2-c]pyridazine"),
    # thieno[3,2-d][1,2,3]triazine (CH at 4,6,7)
    ("c1cc2nnncc2s1",       "thieno[3,2-d][1,2,3]triazine"),
    ("Cc1nnnc2ccsc12",      "4-methylthieno[3,2-d][1,2,3]triazine"),
    ("Cc1cc2nnncc2s1",      "6-methylthieno[3,2-d][1,2,3]triazine"),
    ("Cc1csc2cnnnc12",      "7-methylthieno[3,2-d][1,2,3]triazine"),
    # thieno[3,2-e][1,2,4]triazine (CH at 3,5,6)
    ("c1nnc2sccc2n1",       "thieno[3,2-e][1,2,4]triazine"),
    ("Cc1nnc2sccc2n1",      "3-methylthieno[3,2-e][1,2,4]triazine"),
    ("Cc1csc2nncnc12",      "5-methylthieno[3,2-e][1,2,4]triazine"),
    ("Cc1cc2ncnnc2s1",      "6-methylthieno[3,2-e][1,2,4]triazine"),
    # thieno[3,4-b]pyridine (CH at 2,3,4,5,7)
    ("c1cnc2cscc2c1",       "thieno[3,4-b]pyridine"),
    ("Cc1ccc2cscc2n1",      "2-methylthieno[3,4-b]pyridine"),
    ("Cc1cnc2cscc2c1",      "3-methylthieno[3,4-b]pyridine"),
    ("Cc1ccnc2cscc12",      "4-methylthieno[3,4-b]pyridine"),
    ("Cc1scc2ncccc12",      "5-methylthieno[3,4-b]pyridine"),
    ("Cc1scc2cccnc12",      "7-methylthieno[3,4-b]pyridine"),
    # thieno[3,4-c]pyridine (CH at 1,3,4,6,7)
    ("c1cc2cscc2cn1",       "thieno[3,4-c]pyridine"),
    ("Cc1scc2cnccc12",      "1-methylthieno[3,4-c]pyridine"),
    ("Cc1scc2ccncc12",      "3-methylthieno[3,4-c]pyridine"),
    ("Cc1nccc2cscc12",      "4-methylthieno[3,4-c]pyridine"),
    ("Cc1cc2cscc2cn1",      "6-methylthieno[3,4-c]pyridine"),
    ("Cc1cncc2cscc12",      "7-methylthieno[3,4-c]pyridine"),
    # thieno[3,4-c]pyridazine (CH at 3,4,5,7)
    ("c1cc2cscc2nn1",       "thieno[3,4-c]pyridazine"),
    ("Cc1cc2cscc2nn1",      "3-methylthieno[3,4-c]pyridazine"),
    ("Cc1cnnc2cscc12",      "4-methylthieno[3,4-c]pyridazine"),
    ("Cc1scc2nnccc12",      "5-methylthieno[3,4-c]pyridazine"),
    ("Cc1scc2ccnnc12",      "7-methylthieno[3,4-c]pyridazine"),
    # thieno[3,4-d][1,2,3]triazine (CH at 4,5,7)
    ("c1nnnc2cscc12",       "thieno[3,4-d][1,2,3]triazine"),
    ("Cc1nnnc2cscc12",      "4-methylthieno[3,4-d][1,2,3]triazine"),
    ("Cc1scc2nnncc12",      "5-methylthieno[3,4-d][1,2,3]triazine"),
    ("Cc1scc2cnnnc12",      "7-methylthieno[3,4-d][1,2,3]triazine"),
    # thieno[3,4-d]pyridazine (CH at 1,5)
    ("c1nncc2cscc12",       "thieno[3,4-d]pyridazine"),
    ("Cc1nncc2cscc12",      "1-methylthieno[3,4-d]pyridazine"),
    ("Cc1scc2cnncc12",      "5-methylthieno[3,4-d]pyridazine"),
    # thieno[3,4-e][1,2,4]triazine (CH at 3,5,7)
    ("c1nnc2cscc2n1",       "thieno[3,4-e][1,2,4]triazine"),
    ("Cc1nnc2cscc2n1",      "3-methylthieno[3,4-e][1,2,4]triazine"),
    ("Cc1scc2nncnc12",      "5-methylthieno[3,4-e][1,2,4]triazine"),
    ("Cc1scc2ncnnc12",      "7-methylthieno[3,4-e][1,2,4]triazine"),
    # thieno[3,4-e]pyrazine (CH at 2,5)
    ("c1cnc2cscc2n1",       "thieno[3,4-e]pyrazine"),
    ("Cc1cnc2cscc2n1",      "2-methylthieno[3,4-e]pyrazine"),
    ("Cc1scc2nccnc12",      "5-methylthieno[3,4-e]pyrazine"),
    # thieno[2,3-b]thiophene (CH at 2,3)
    ("c1cc2ccsc2s1",        "thieno[2,3-b]thiophene"),
    ("Cc1cc2ccsc2s1",       "2-methylthieno[2,3-b]thiophene"),
    ("Cc1csc2sccc12",       "3-methylthieno[2,3-b]thiophene"),
    # thieno[3,2-b]thiophene (CH at 2,3)
    ("c1cc2sccc2s1",        "thieno[3,2-b]thiophene"),
    ("Cc1cc2sccc2s1",       "2-methylthieno[3,2-b]thiophene"),
    ("Cc1csc2ccsc12",       "3-methylthieno[3,2-b]thiophene"),
])
def test_phase714(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
