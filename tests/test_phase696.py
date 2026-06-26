"""Phase 696: methyl derivatives of c-fused 9-atom bicyclics with pyridazine
(1H-imidazo[4,5-c]pyridazine, 1H-pyrazolo[3,4-c]pyridazine,
1H-pyrazolo[4,5-c]pyridazine, 1H/2H-[1,2,3]triazolo[4,5-c]pyridazine,
1H-[1,2,3]triazolo[5,4-c]pyridazine).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-imidazo[4,5-c]pyridazine (CH at 3,4,6)
    ("c1cc2[nH]cnc2nn1",   "1H-imidazo[4,5-c]pyridazine"),
    ("Cc1cc2[nH]cnc2nn1",  "3-methyl-1H-imidazo[4,5-c]pyridazine"),
    ("Cc1cnnc2nc[nH]c12",  "4-methyl-1H-imidazo[4,5-c]pyridazine"),
    ("Cc1nc2nnccc2[nH]1",  "6-methyl-1H-imidazo[4,5-c]pyridazine"),
    # 1H-pyrazolo[3,4-c]pyridazine (CH at 3,4,5)
    ("c1cc2c[nH]nc2nn1",   "1H-pyrazolo[3,4-c]pyridazine"),
    ("Cc1[nH]nc2nnccc12",  "3-methyl-1H-pyrazolo[3,4-c]pyridazine"),
    ("Cc1cnnc2n[nH]cc12",  "4-methyl-1H-pyrazolo[3,4-c]pyridazine"),
    ("Cc1cc2c[nH]nc2nn1",  "5-methyl-1H-pyrazolo[3,4-c]pyridazine"),
    # 1H-pyrazolo[4,5-c]pyridazine (CH at 3,6,7)
    ("c1cc2[nH]ncc2nn1",   "1H-pyrazolo[4,5-c]pyridazine"),
    ("Cc1n[nH]c2ccnnc12",  "3-methyl-1H-pyrazolo[4,5-c]pyridazine"),
    ("Cc1cc2[nH]ncc2nn1",  "6-methyl-1H-pyrazolo[4,5-c]pyridazine"),
    ("Cc1cnnc2cn[nH]c12",  "7-methyl-1H-pyrazolo[4,5-c]pyridazine"),
    # 1H-[1,2,3]triazolo[4,5-c]pyridazine (CH at 6,7)
    ("c1cc2[nH]nnc2nn1",   "1H-[1,2,3]triazolo[4,5-c]pyridazine"),
    ("Cc1cc2[nH]nnc2nn1",  "6-methyl-1H-[1,2,3]triazolo[4,5-c]pyridazine"),
    ("Cc1cnnc2nn[nH]c12",  "7-methyl-1H-[1,2,3]triazolo[4,5-c]pyridazine"),
    # 2H-[1,2,3]triazolo[4,5-c]pyridazine (CH at 6,7)
    ("c1cc2n[nH]nc2nn1",   "2H-[1,2,3]triazolo[4,5-c]pyridazine"),
    ("Cc1cc2n[nH]nc2nn1",  "6-methyl-2H-[1,2,3]triazolo[4,5-c]pyridazine"),
    ("Cc1cnnc2n[nH]nc12",  "7-methyl-2H-[1,2,3]triazolo[4,5-c]pyridazine"),
    # 1H-[1,2,3]triazolo[5,4-c]pyridazine (CH at 6,7)
    ("c1cc2nn[nH]c2nn1",   "1H-[1,2,3]triazolo[5,4-c]pyridazine"),
    ("Cc1cc2nn[nH]c2nn1",  "6-methyl-1H-[1,2,3]triazolo[5,4-c]pyridazine"),
    ("Cc1cnnc2[nH]nnc12",  "7-methyl-1H-[1,2,3]triazolo[5,4-c]pyridazine"),
])
def test_phase696(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
