"""Phase 698: methyl derivatives of aromatic 1H-pyrrolo[x,y-d/e] bicyclics
(pyrimidine, pyridazine, and pyrazine partners; 4 parent compounds).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 7H-pyrrolo[2,3-d]pyrimidine (NH at N7; Phase 838 correction)
    ("c1ncc2cc[nH]c2n1",    "7H-pyrrolo[2,3-d]pyrimidine"),
    ("Cc1ncc2cc[nH]c2n1",   "2-methyl-7H-pyrrolo[2,3-d]pyrimidine"),
    ("Cc1ncnc2[nH]ccc12",   "4-methyl-7H-pyrrolo[2,3-d]pyrimidine"),
    ("Cc1c[nH]c2ncncc12",   "5-methyl-7H-pyrrolo[2,3-d]pyrimidine"),
    ("Cc1cc2cncnc2[nH]1",   "6-methyl-7H-pyrrolo[2,3-d]pyrimidine"),
    # 1H-pyrrolo[3,2-d]pyrimidine (CH at 2,4,6,7)
    ("c1cc2[nH]cncc-2n1",   "1H-pyrrolo[3,2-d]pyrimidine"),
    ("Cc1ncc2nccc-2[nH]1",  "2-methyl-1H-pyrrolo[3,2-d]pyrimidine"),
    ("Cc1nc[nH]c2ccnc1-2",  "4-methyl-1H-pyrrolo[3,2-d]pyrimidine"),
    ("Cc1cc2[nH]cncc-2n1",  "6-methyl-1H-pyrrolo[3,2-d]pyrimidine"),
    ("Cc1cnc2cnc[nH]c1-2",  "7-methyl-1H-pyrrolo[3,2-d]pyrimidine"),
    # 1H-pyrrolo[2,3-d]pyridazine (CH at 2,3,4,7)
    ("c1cc2cnncc2[nH]1",    "1H-pyrrolo[2,3-d]pyridazine"),
    ("Cc1cc2cnncc2[nH]1",   "2-methyl-1H-pyrrolo[2,3-d]pyridazine"),
    ("Cc1c[nH]c2cnncc12",   "3-methyl-1H-pyrrolo[2,3-d]pyridazine"),
    ("Cc1nncc2[nH]ccc12",   "4-methyl-1H-pyrrolo[2,3-d]pyridazine"),
    ("Cc1nncc2cc[nH]c12",   "7-methyl-1H-pyrrolo[2,3-d]pyridazine"),
    # 1H-pyrrolo[2,3-e]pyrazine (CH at 2,3,6,7)
    ("c1c[nH]c2ccnc-2n1",   "1H-pyrrolo[2,3-e]pyrazine"),
    ("Cc1cnc2nccc-2[nH]1",  "2-methyl-1H-pyrrolo[2,3-e]pyrazine"),
    ("Cc1c[nH]c2ccnc-2n1",  "3-methyl-1H-pyrrolo[2,3-e]pyrazine"),
    ("Cc1cc2[nH]ccnc-2n1",  "6-methyl-1H-pyrrolo[2,3-e]pyrazine"),
    ("Cc1cnc2ncc[nH]c1-2",  "7-methyl-1H-pyrrolo[2,3-e]pyrazine"),
])
def test_phase698(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
