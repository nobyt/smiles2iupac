"""Phase 697: methyl derivatives of 1H-pyrazolo[x,y-d/e] bicyclics
(pyrimidine, pyridazine, and pyrazine partners; 8 parent compounds).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-pyrazolo[3,4-d]pyrimidine (CH at 3,4,6)
    ("c1ncc2c[nH]nc2n1",   "1H-pyrazolo[3,4-d]pyrimidine"),
    ("Cc1[nH]nc2ncncc12",  "3-methyl-1H-pyrazolo[3,4-d]pyrimidine"),
    ("Cc1ncnc2n[nH]cc12",  "4-methyl-1H-pyrazolo[3,4-d]pyrimidine"),
    ("Cc1ncc2c[nH]nc2n1",  "6-methyl-1H-pyrazolo[3,4-d]pyrimidine"),
    # 1H-pyrazolo[4,5-d]pyrimidine (CH at 3,5,7)
    ("c1ncc2[nH]ncc2n1",   "1H-pyrazolo[4,5-d]pyrimidine"),
    ("Cc1n[nH]c2cncnc12",  "3-methyl-1H-pyrazolo[4,5-d]pyrimidine"),
    ("Cc1ncc2[nH]ncc2n1",  "5-methyl-1H-pyrazolo[4,5-d]pyrimidine"),
    ("Cc1ncnc2cn[nH]c12",  "7-methyl-1H-pyrazolo[4,5-d]pyrimidine"),
    # 1H-pyrazolo[4,3-d]pyrimidine (CH at 3,5,7)
    ("c1ncc2n[nH]cc2n1",   "1H-pyrazolo[4,3-d]pyrimidine"),
    ("Cc1[nH]nc2cncnc12",  "3-methyl-1H-pyrazolo[4,3-d]pyrimidine"),
    ("Cc1ncc2n[nH]cc2n1",  "5-methyl-1H-pyrazolo[4,3-d]pyrimidine"),
    ("Cc1ncnc2c[nH]nc12",  "7-methyl-1H-pyrazolo[4,3-d]pyrimidine"),
    # 1H-pyrazolo[5,4-d]pyrimidine (CH at 3,4,6)
    ("c1ncc2cn[nH]c2n1",   "1H-pyrazolo[5,4-d]pyrimidine"),
    ("Cc1n[nH]c2ncncc12",  "3-methyl-1H-pyrazolo[5,4-d]pyrimidine"),
    ("Cc1ncnc2[nH]ncc12",  "4-methyl-1H-pyrazolo[5,4-d]pyrimidine"),
    ("Cc1ncc2cn[nH]c2n1",  "6-methyl-1H-pyrazolo[5,4-d]pyrimidine"),
    # 1H-pyrazolo[3,4-d]pyridazine (CH at 3,4,7)
    ("c1nncc2n[nH]cc12",   "1H-pyrazolo[3,4-d]pyridazine"),
    ("Cc1[nH]nc2cnncc12",  "3-methyl-1H-pyrazolo[3,4-d]pyridazine"),
    ("Cc1nncc2n[nH]cc12",  "4-methyl-1H-pyrazolo[3,4-d]pyridazine"),
    ("Cc1nncc2c[nH]nc12",  "7-methyl-1H-pyrazolo[3,4-d]pyridazine"),
    # 1H-pyrazolo[4,5-d]pyridazine (CH at 3,4,7)
    ("c1nncc2[nH]ncc12",   "1H-pyrazolo[4,5-d]pyridazine"),
    ("Cc1n[nH]c2cnncc12",  "3-methyl-1H-pyrazolo[4,5-d]pyridazine"),
    ("Cc1nncc2[nH]ncc12",  "4-methyl-1H-pyrazolo[4,5-d]pyridazine"),
    ("Cc1nncc2cn[nH]c12",  "7-methyl-1H-pyrazolo[4,5-d]pyridazine"),
    # 1H-pyrazolo[3,4-e]pyrazine (CH at 3,5,6)
    ("c1cnc2n[nH]cc2n1",   "1H-pyrazolo[3,4-e]pyrazine"),
    ("Cc1[nH]nc2nccnc12",  "3-methyl-1H-pyrazolo[3,4-e]pyrazine"),
    ("Cc1cnc2n[nH]cc2n1",  "5-methyl-1H-pyrazolo[3,4-e]pyrazine"),
    ("Cc1cnc2c[nH]nc2n1",  "6-methyl-1H-pyrazolo[3,4-e]pyrazine"),
    # 1H-pyrazolo[4,5-e]pyrazine (CH at 3,5,6)
    ("c1cnc2[nH]ncc2n1",   "1H-pyrazolo[4,5-e]pyrazine"),
    ("Cc1n[nH]c2nccnc12",  "3-methyl-1H-pyrazolo[4,5-e]pyrazine"),
    ("Cc1cnc2[nH]ncc2n1",  "5-methyl-1H-pyrazolo[4,5-e]pyrazine"),
    ("Cc1cnc2cn[nH]c2n1",  "6-methyl-1H-pyrazolo[4,5-e]pyrazine"),
])
def test_phase697(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
