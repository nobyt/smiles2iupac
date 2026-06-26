"""Phase 709: methyl derivatives of isothiazolo[x,y-b/c/d/e] series
(28 parent compounds: all b/c/d/e ring-fusion isomers with pyridine,
pyridazine, pyrimidine, pyrazine, and triazine partners).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isothiazolo[3,4-b]pyridine (CH at 3,4,5,6)
    ("c1cnc2nscc2c1",       "isothiazolo[3,4-b]pyridine"),
    ("Cc1snc2ncccc12",      "3-methylisothiazolo[3,4-b]pyridine"),
    ("Cc1ccnc2nscc12",      "4-methylisothiazolo[3,4-b]pyridine"),
    ("Cc1cnc2nscc2c1",      "5-methylisothiazolo[3,4-b]pyridine"),
    ("Cc1ccc2csnc2n1",      "6-methylisothiazolo[3,4-b]pyridine"),
    # isothiazolo[3,4-c]pyridazine (CH at 3,4,5)
    ("c1cc2csnc2nn1",       "isothiazolo[3,4-c]pyridazine"),
    ("Cc1snc2nnccc12",      "3-methylisothiazolo[3,4-c]pyridazine"),
    ("Cc1cnnc2nscc12",      "4-methylisothiazolo[3,4-c]pyridazine"),
    ("Cc1cc2csnc2nn1",      "5-methylisothiazolo[3,4-c]pyridazine"),
    # isothiazolo[3,4-c]pyridine (CH at 3,4,5,7)
    ("c1cc2csnc2cn1",       "isothiazolo[3,4-c]pyridine"),
    ("Cc1snc2cnccc12",      "3-methylisothiazolo[3,4-c]pyridine"),
    ("Cc1cncc2nscc12",      "4-methylisothiazolo[3,4-c]pyridine"),
    ("Cc1cc2csnc2cn1",      "5-methylisothiazolo[3,4-c]pyridine"),
    ("Cc1nccc2csnc12",      "7-methylisothiazolo[3,4-c]pyridine"),
    # isothiazolo[3,4-d][1,2,3]triazine (CH at 4,5)
    ("c1nnnc2nscc12",       "isothiazolo[3,4-d][1,2,3]triazine"),
    ("Cc1nnnc2nscc12",      "4-methylisothiazolo[3,4-d][1,2,3]triazine"),
    ("Cc1snc2nnncc12",      "5-methylisothiazolo[3,4-d][1,2,3]triazine"),
    # isothiazolo[3,4-d]pyridazine (CH at 3,4,7)
    ("c1nncc2nscc12",       "isothiazolo[3,4-d]pyridazine"),
    ("Cc1snc2cnncc12",      "3-methylisothiazolo[3,4-d]pyridazine"),
    ("Cc1nncc2nscc12",      "4-methylisothiazolo[3,4-d]pyridazine"),
    ("Cc1nncc2csnc12",      "7-methylisothiazolo[3,4-d]pyridazine"),
    # isothiazolo[3,4-d]pyrimidine (CH at 3,4,6)
    ("c1ncc2csnc2n1",       "isothiazolo[3,4-d]pyrimidine"),
    ("Cc1snc2ncncc12",      "3-methylisothiazolo[3,4-d]pyrimidine"),
    ("Cc1ncnc2nscc12",      "4-methylisothiazolo[3,4-d]pyrimidine"),
    ("Cc1ncc2csnc2n1",      "6-methylisothiazolo[3,4-d]pyrimidine"),
    # isothiazolo[3,4-e][1,2,4]triazine (CH at 3,7)
    ("c1nnc2csnc2n1",       "isothiazolo[3,4-e][1,2,4]triazine"),
    ("Cc1nnc2csnc2n1",      "3-methylisothiazolo[3,4-e][1,2,4]triazine"),
    ("Cc1snc2ncnnc12",      "7-methylisothiazolo[3,4-e][1,2,4]triazine"),
    # isothiazolo[3,4-e]pyrazine (CH at 3,5,6)
    ("c1cnc2nscc2n1",       "isothiazolo[3,4-e]pyrazine"),
    ("Cc1snc2nccnc12",      "3-methylisothiazolo[3,4-e]pyrazine"),
    ("Cc1cnc2nscc2n1",      "5-methylisothiazolo[3,4-e]pyrazine"),
    ("Cc1cnc2csnc2n1",      "6-methylisothiazolo[3,4-e]pyrazine"),
    # isothiazolo[4,3-b]pyridine (CH at 3,5,6,7)
    ("c1cnc2csnc2c1",       "isothiazolo[4,3-b]pyridine"),
    ("Cc1snc2cccnc12",      "3-methylisothiazolo[4,3-b]pyridine"),
    ("Cc1ccc2nscc2n1",      "5-methylisothiazolo[4,3-b]pyridine"),
    ("Cc1cnc2csnc2c1",      "6-methylisothiazolo[4,3-b]pyridine"),
    ("Cc1ccnc2csnc12",      "7-methylisothiazolo[4,3-b]pyridine"),
    # isothiazolo[4,3-c]pyridazine (CH at 3,6,7)
    ("c1cc2nscc2nn1",       "isothiazolo[4,3-c]pyridazine"),
    ("Cc1snc2ccnnc12",      "3-methylisothiazolo[4,3-c]pyridazine"),
    ("Cc1cc2nscc2nn1",      "6-methylisothiazolo[4,3-c]pyridazine"),
    ("Cc1cnnc2csnc12",      "7-methylisothiazolo[4,3-c]pyridazine"),
    # isothiazolo[4,3-c]pyridine (CH at 3,4,6,7)
    ("c1cc2nscc2cn1",       "isothiazolo[4,3-c]pyridine"),
    ("Cc1snc2ccncc12",      "3-methylisothiazolo[4,3-c]pyridine"),
    ("Cc1nccc2nscc12",      "4-methylisothiazolo[4,3-c]pyridine"),
    ("Cc1cc2nscc2cn1",      "6-methylisothiazolo[4,3-c]pyridine"),
    ("Cc1cncc2csnc12",      "7-methylisothiazolo[4,3-c]pyridine"),
    # isothiazolo[4,3-d][1,2,3]triazine (CH at 4,7)
    ("c1snc2cnnnc12",       "isothiazolo[4,3-d][1,2,3]triazine"),
    ("Cc1nnnc2csnc12",      "4-methylisothiazolo[4,3-d][1,2,3]triazine"),
    ("Cc1snc2cnnnc12",      "7-methylisothiazolo[4,3-d][1,2,3]triazine"),
    # isothiazolo[4,3-d]pyrimidine (CH at 3,5,7)
    ("c1ncc2nscc2n1",       "isothiazolo[4,3-d]pyrimidine"),
    ("Cc1snc2cncnc12",      "3-methylisothiazolo[4,3-d]pyrimidine"),
    ("Cc1ncc2nscc2n1",      "5-methylisothiazolo[4,3-d]pyrimidine"),
    ("Cc1ncnc2csnc12",      "7-methylisothiazolo[4,3-d]pyrimidine"),
    # isothiazolo[4,3-e][1,2,4]triazine (CH at 3,5)
    ("c1nnc2nscc2n1",       "isothiazolo[4,3-e][1,2,4]triazine"),
    ("Cc1snc2nncnc12",      "3-methylisothiazolo[4,3-e][1,2,4]triazine"),
    ("Cc1nnc2nscc2n1",      "5-methylisothiazolo[4,3-e][1,2,4]triazine"),
    # isothiazolo[4,5-b]pyridine (CH at 3,5,6,7)
    ("c1cnc2cnsc2c1",       "isothiazolo[4,5-b]pyridine"),
    ("Cc1nsc2cccnc12",      "3-methylisothiazolo[4,5-b]pyridine"),
    ("Cc1ccc2sncc2n1",      "5-methylisothiazolo[4,5-b]pyridine"),
    ("Cc1cnc2cnsc2c1",      "6-methylisothiazolo[4,5-b]pyridine"),
    ("Cc1ccnc2cnsc12",      "7-methylisothiazolo[4,5-b]pyridine"),
    # isothiazolo[4,5-c]pyridazine (CH at 3,6,7)
    ("c1cc2sncc2nn1",       "isothiazolo[4,5-c]pyridazine"),
    ("Cc1nsc2ccnnc12",      "3-methylisothiazolo[4,5-c]pyridazine"),
    ("Cc1cc2sncc2nn1",      "6-methylisothiazolo[4,5-c]pyridazine"),
    ("Cc1cnnc2cnsc12",      "7-methylisothiazolo[4,5-c]pyridazine"),
    # isothiazolo[4,5-c]pyridine (CH at 3,4,6,7)
    ("c1cc2sncc2cn1",       "isothiazolo[4,5-c]pyridine"),
    ("Cc1nsc2ccncc12",      "3-methylisothiazolo[4,5-c]pyridine"),
    ("Cc1nccc2sncc12",      "4-methylisothiazolo[4,5-c]pyridine"),
    ("Cc1cc2sncc2cn1",      "6-methylisothiazolo[4,5-c]pyridine"),
    ("Cc1cncc2cnsc12",      "7-methylisothiazolo[4,5-c]pyridine"),
    # isothiazolo[4,5-d][1,2,3]triazine (CH at 4,7)
    ("c1nsc2cnnnc12",       "isothiazolo[4,5-d][1,2,3]triazine"),
    ("Cc1nnnc2cnsc12",      "4-methylisothiazolo[4,5-d][1,2,3]triazine"),
    ("Cc1nsc2cnnnc12",      "7-methylisothiazolo[4,5-d][1,2,3]triazine"),
    # isothiazolo[4,5-d]pyrimidine (CH at 3,5,7)
    ("c1ncc2sncc2n1",       "isothiazolo[4,5-d]pyrimidine"),
    ("Cc1nsc2cncnc12",      "3-methylisothiazolo[4,5-d]pyrimidine"),
    ("Cc1ncc2sncc2n1",      "5-methylisothiazolo[4,5-d]pyrimidine"),
    ("Cc1ncnc2cnsc12",      "7-methylisothiazolo[4,5-d]pyrimidine"),
    # isothiazolo[4,5-e][1,2,4]triazine (CH at 3,5)
    ("c1nnc2sncc2n1",       "isothiazolo[4,5-e][1,2,4]triazine"),
    ("Cc1nsc2nncnc12",      "3-methylisothiazolo[4,5-e][1,2,4]triazine"),
    ("Cc1nnc2sncc2n1",      "5-methylisothiazolo[4,5-e][1,2,4]triazine"),
    # isothiazolo[4,5-e]pyrazine (CH at 3,5,6)
    ("c1cnc2sncc2n1",       "isothiazolo[4,5-e]pyrazine"),
    ("Cc1nsc2nccnc12",      "3-methylisothiazolo[4,5-e]pyrazine"),
    ("Cc1cnc2sncc2n1",      "5-methylisothiazolo[4,5-e]pyrazine"),
    ("Cc1cnc2cnsc2n1",      "6-methylisothiazolo[4,5-e]pyrazine"),
    # isothiazolo[5,4-b]pyridine (CH at 3,4,5,6)
    ("c1cnc2sncc2c1",       "isothiazolo[5,4-b]pyridine"),
    ("Cc1nsc2ncccc12",      "3-methylisothiazolo[5,4-b]pyridine"),
    ("Cc1ccnc2sncc12",      "4-methylisothiazolo[5,4-b]pyridine"),
    ("Cc1cnc2sncc2c1",      "5-methylisothiazolo[5,4-b]pyridine"),
    ("Cc1ccc2cnsc2n1",      "6-methylisothiazolo[5,4-b]pyridine"),
    # isothiazolo[5,4-c]pyridazine (CH at 3,4,5)
    ("c1cc2cnsc2nn1",       "isothiazolo[5,4-c]pyridazine"),
    ("Cc1nsc2nnccc12",      "3-methylisothiazolo[5,4-c]pyridazine"),
    ("Cc1cnnc2sncc12",      "4-methylisothiazolo[5,4-c]pyridazine"),
    ("Cc1cc2cnsc2nn1",      "5-methylisothiazolo[5,4-c]pyridazine"),
    # isothiazolo[5,4-c]pyridine (CH at 3,4,5,7)
    ("c1cc2cnsc2cn1",       "isothiazolo[5,4-c]pyridine"),
    ("Cc1nsc2cnccc12",      "3-methylisothiazolo[5,4-c]pyridine"),
    ("Cc1cncc2sncc12",      "4-methylisothiazolo[5,4-c]pyridine"),
    ("Cc1cc2cnsc2cn1",      "5-methylisothiazolo[5,4-c]pyridine"),
    ("Cc1nccc2cnsc12",      "7-methylisothiazolo[5,4-c]pyridine"),
    # isothiazolo[5,4-d][1,2,3]triazine (CH at 4,5)
    ("c1nnnc2sncc12",       "isothiazolo[5,4-d][1,2,3]triazine"),
    ("Cc1nnnc2sncc12",      "4-methylisothiazolo[5,4-d][1,2,3]triazine"),
    ("Cc1nsc2nnncc12",      "5-methylisothiazolo[5,4-d][1,2,3]triazine"),
    # isothiazolo[5,4-d]pyridazine (CH at 3,4,7)
    ("c1nncc2sncc12",       "isothiazolo[5,4-d]pyridazine"),
    ("Cc1nsc2cnncc12",      "3-methylisothiazolo[5,4-d]pyridazine"),
    ("Cc1nncc2sncc12",      "4-methylisothiazolo[5,4-d]pyridazine"),
    ("Cc1nncc2cnsc12",      "7-methylisothiazolo[5,4-d]pyridazine"),
    # isothiazolo[5,4-d]pyrimidine (CH at 3,4,6)
    ("c1ncc2cnsc2n1",       "isothiazolo[5,4-d]pyrimidine"),
    ("Cc1nsc2ncncc12",      "3-methylisothiazolo[5,4-d]pyrimidine"),
    ("Cc1ncnc2sncc12",      "4-methylisothiazolo[5,4-d]pyrimidine"),
    ("Cc1ncc2cnsc2n1",      "6-methylisothiazolo[5,4-d]pyrimidine"),
    # isothiazolo[5,4-e][1,2,4]triazine (CH at 3,7)
    ("c1nnc2cnsc2n1",       "isothiazolo[5,4-e][1,2,4]triazine"),
    ("Cc1nnc2cnsc2n1",      "3-methylisothiazolo[5,4-e][1,2,4]triazine"),
    ("Cc1nsc2ncnnc12",      "7-methylisothiazolo[5,4-e][1,2,4]triazine"),
])
def test_phase709(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
