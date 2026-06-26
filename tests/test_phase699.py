"""Phase 699: methyl derivatives of 1H-pyrrolo[x,y-d/e] with triazine partners
and 1H-pyrrolo[3,4-d/e] series (8 parent compounds).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-pyrrolo[2,3-d][1,2,3]triazine (CH at 4,5,6)
    ("c1cc2cnnnc2[nH]1",    "1H-pyrrolo[2,3-d][1,2,3]triazine"),
    ("Cc1nnnc2[nH]ccc12",   "4-methyl-1H-pyrrolo[2,3-d][1,2,3]triazine"),
    ("Cc1c[nH]c2nnncc12",   "5-methyl-1H-pyrrolo[2,3-d][1,2,3]triazine"),
    ("Cc1cc2cnnnc2[nH]1",   "6-methyl-1H-pyrrolo[2,3-d][1,2,3]triazine"),
    # 1H-pyrrolo[3,2-d][1,2,3]triazine (CH at 4,6,7)
    ("c1cc2[nH]nncc-2n1",   "1H-pyrrolo[3,2-d][1,2,3]triazine"),
    ("Cc1nn[nH]c2ccnc1-2",  "4-methyl-1H-pyrrolo[3,2-d][1,2,3]triazine"),
    ("Cc1cc2[nH]nncc-2n1",  "6-methyl-1H-pyrrolo[3,2-d][1,2,3]triazine"),
    ("Cc1cnc2cnn[nH]c1-2",  "7-methyl-1H-pyrrolo[3,2-d][1,2,3]triazine"),
    # 1H-pyrrolo[2,3-e][1,2,4]triazine (CH at 3,6,7)
    ("c1cc2[nH]ncnc-2n1",   "1H-pyrrolo[2,3-e][1,2,4]triazine"),
    ("Cc1n[nH]c2ccnc-2n1",  "3-methyl-1H-pyrrolo[2,3-e][1,2,4]triazine"),
    ("Cc1cc2[nH]ncnc-2n1",  "6-methyl-1H-pyrrolo[2,3-e][1,2,4]triazine"),
    ("Cc1cnc2ncn[nH]c1-2",  "7-methyl-1H-pyrrolo[2,3-e][1,2,4]triazine"),
    # 1H-pyrrolo[3,2-e][1,2,4]triazine (CH at 3,5,6)
    ("c1nnc2[nH]ccc2n1",    "1H-pyrrolo[3,2-e][1,2,4]triazine"),
    ("Cc1nnc2[nH]ccc2n1",   "3-methyl-1H-pyrrolo[3,2-e][1,2,4]triazine"),
    ("Cc1c[nH]c2nncnc12",   "5-methyl-1H-pyrrolo[3,2-e][1,2,4]triazine"),
    ("Cc1cc2ncnnc2[nH]1",   "6-methyl-1H-pyrrolo[3,2-e][1,2,4]triazine"),
    # 1H-pyrrolo[3,4-d]pyrimidine (CH at 2,4,5,7)
    ("c1ncc2cncc-2[nH]1",   "1H-pyrrolo[3,4-d]pyrimidine"),
    ("Cc1ncc2cncc-2[nH]1",  "2-methyl-1H-pyrrolo[3,4-d]pyrimidine"),
    ("Cc1nc[nH]c2cncc1-2",  "4-methyl-1H-pyrrolo[3,4-d]pyrimidine"),
    ("Cc1ncc2[nH]cncc1-2",  "5-methyl-1H-pyrrolo[3,4-d]pyrimidine"),
    ("Cc1ncc2cnc[nH]c1-2",  "7-methyl-1H-pyrrolo[3,4-d]pyrimidine"),
    # 1H-pyrrolo[3,4-e]pyrazine (CH at 2,3,5,7)
    ("c1c[nH]c2cncc-2n1",   "1H-pyrrolo[3,4-e]pyrazine"),
    ("Cc1cnc2cncc-2[nH]1",  "2-methyl-1H-pyrrolo[3,4-e]pyrazine"),
    ("Cc1c[nH]c2cncc-2n1",  "3-methyl-1H-pyrrolo[3,4-e]pyrazine"),
    ("Cc1ncc2[nH]ccnc1-2",  "5-methyl-1H-pyrrolo[3,4-e]pyrazine"),
    ("Cc1ncc2ncc[nH]c1-2",  "7-methyl-1H-pyrrolo[3,4-e]pyrazine"),
    # 1H-pyrrolo[3,4-d][1,2,3]triazine (CH at 4,5,7)
    ("c1ncc2[nH]nncc1-2",   "1H-pyrrolo[3,4-d][1,2,3]triazine"),
    ("Cc1nn[nH]c2cncc1-2",  "4-methyl-1H-pyrrolo[3,4-d][1,2,3]triazine"),
    ("Cc1ncc2[nH]nncc1-2",  "5-methyl-1H-pyrrolo[3,4-d][1,2,3]triazine"),
    ("Cc1ncc2cnn[nH]c1-2",  "7-methyl-1H-pyrrolo[3,4-d][1,2,3]triazine"),
    # 1H-pyrrolo[3,4-e][1,2,4]triazine (CH at 3,5,7)
    ("c1n[nH]c2cncc-2n1",   "1H-pyrrolo[3,4-e][1,2,4]triazine"),
    ("Cc1n[nH]c2cncc-2n1",  "3-methyl-1H-pyrrolo[3,4-e][1,2,4]triazine"),
    ("Cc1ncc2[nH]ncnc1-2",  "5-methyl-1H-pyrrolo[3,4-e][1,2,4]triazine"),
    ("Cc1ncc2ncn[nH]c1-2",  "7-methyl-1H-pyrrolo[3,4-e][1,2,4]triazine"),
])
def test_phase699(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
