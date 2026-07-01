"""Phase 729: missing methyl derivatives of imidazo/triazolo-fused pyridines,
pyrazines, pyridazines, pyrimidines, and triazines; plus monocyclic
1,2,4-triazine, 1,2-thiazole (isothiazole PIN), 1,3-oxazole (oxazole PIN), and 1,3-thiazole (thiazole PIN).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-imidazo[4,5-b]pyridine (CH at 2,5,6,7; two tautomers)
    ("c1cnc2[nH]cnc2c1",        "1H-imidazo[4,5-b]pyridine"),
    ("Cc1nc2cccnc2[nH]1",       "2-methyl-1H-imidazo[4,5-b]pyridine"),
    ("Cc1ccc2nc[nH]c2n1",       "5-methyl-1H-imidazo[4,5-b]pyridine"),
    ("Cc1cnc2[nH]cnc2c1",       "6-methyl-1H-imidazo[4,5-b]pyridine"),
    ("Cc1ccnc2[nH]cnc12",       "7-methyl-1H-imidazo[4,5-b]pyridine"),
    # 3H-imidazo[4,5-b]pyridine (tautomer; CH at 2,5,6,7)
    ("c1cnc2nc[nH]c2c1",        "3H-imidazo[4,5-b]pyridine"),
    ("Cc1nc2ncccc2[nH]1",       "2-methyl-3H-imidazo[4,5-b]pyridine"),
    ("Cc1ccc2[nH]cnc2n1",       "5-methyl-3H-imidazo[4,5-b]pyridine"),
    ("Cc1cnc2nc[nH]c2c1",       "6-methyl-3H-imidazo[4,5-b]pyridine"),
    ("Cc1ccnc2nc[nH]c12",       "7-methyl-3H-imidazo[4,5-b]pyridine"),
    # 1H-imidazo[4,5-c]pyridine (CH at 2,4,6,7)
    ("c1cc2[nH]cnc2cn1",        "1H-imidazo[4,5-c]pyridine"),
    ("Cc1nc2cnccc2[nH]1",       "2-methyl-1H-imidazo[4,5-c]pyridine"),
    ("Cc1nccc2[nH]cnc12",       "4-methyl-1H-imidazo[4,5-c]pyridine"),
    ("Cc1cc2[nH]cnc2cn1",       "6-methyl-1H-imidazo[4,5-c]pyridine"),
    ("Cc1cncc2nc[nH]c12",       "7-methyl-1H-imidazo[4,5-c]pyridine"),
    # 1H-imidazo[4,5-e]pyrazine (CH at 2,5)
    ("c1cnc2[nH]cnc2n1",        "1H-imidazo[4,5-e]pyrazine"),
    ("Cc1nc2nccnc2[nH]1",       "2-methyl-1H-imidazo[4,5-e]pyrazine"),
    ("Cc1cnc2[nH]cnc2n1",       "5-methyl-1H-imidazo[4,5-e]pyrazine"),
    # 1H-imidazo[4,5-d]pyridazine (CH at 2,4)
    ("c1nc2cnncc2[nH]1",        "1H-imidazo[4,5-d]pyridazine"),
    ("Cc1nc2cnncc2[nH]1",       "2-methyl-1H-imidazo[4,5-d]pyridazine"),
    ("Cc1nncc2[nH]cnc12",       "4-methyl-1H-imidazo[4,5-d]pyridazine"),
    # 1H-imidazo[4,5-d][1,2,3]triazine (CH at 4,6)
    ("c1nc2nnncc2[nH]1",        "1H-imidazo[4,5-d][1,2,3]triazine"),
    ("Cc1nnnc2nc[nH]c12",       "4-methyl-1H-imidazo[4,5-d][1,2,3]triazine"),
    ("Cc1nc2nnncc2[nH]1",       "6-methyl-1H-imidazo[4,5-d][1,2,3]triazine"),
    # 1H-imidazo[5,4-d][1,2,3]triazine (CH at 4,6)
    ("c1nc2cnnnc2[nH]1",        "1H-imidazo[5,4-d][1,2,3]triazine"),
    ("Cc1nnnc2[nH]cnc12",       "4-methyl-1H-imidazo[5,4-d][1,2,3]triazine"),
    ("Cc1nc2cnnnc2[nH]1",       "6-methyl-1H-imidazo[5,4-d][1,2,3]triazine"),
    # 1H-imidazo[4,5-e][1,2,4]triazine (CH at 3,6)
    ("c1nnc2[nH]cnc2n1",        "1H-imidazo[4,5-e][1,2,4]triazine"),
    ("Cc1nnc2[nH]cnc2n1",       "3-methyl-1H-imidazo[4,5-e][1,2,4]triazine"),
    ("Cc1nc2ncnnc2[nH]1",       "6-methyl-1H-imidazo[4,5-e][1,2,4]triazine"),
    # 1H-imidazo[5,4-e][1,2,4]triazine (CH at 3,6)
    ("c1nnc2nc[nH]c2n1",        "1H-imidazo[5,4-e][1,2,4]triazine"),
    ("Cc1nnc2nc[nH]c2n1",       "3-methyl-1H-imidazo[5,4-e][1,2,4]triazine"),
    ("Cc1nc2nncnc2[nH]1",       "6-methyl-1H-imidazo[5,4-e][1,2,4]triazine"),
    # 1H-[1,2,3]triazolo[4,5-b]pyridine (CH at 5,6,7; two tautomers)
    ("c1cnc2[nH]nnc2c1",        "1H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("Cc1ccc2nn[nH]c2n1",       "5-methyl-1H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("Cc1cnc2[nH]nnc2c1",       "6-methyl-1H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("Cc1ccnc2[nH]nnc12",       "7-methyl-1H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("c1cnc2nn[nH]c2c1",        "1H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("Cc1ccc2[nH]nnc2n1",       "5-methyl-1H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("Cc1cnc2nn[nH]c2c1",       "6-methyl-1H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("Cc1ccnc2nn[nH]c12",       "7-methyl-1H-[1,2,3]triazolo[4,5-b]pyridine"),
    # 2H-[1,2,3]triazolo[4,5-b]pyridine (CH at 5,6,7)
    ("c1cnc2n[nH]nc2c1",        "2H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("Cc1ccc2n[nH]nc2n1",       "5-methyl-2H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("Cc1cnc2n[nH]nc2c1",       "6-methyl-2H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("Cc1ccnc2n[nH]nc12",       "7-methyl-2H-[1,2,3]triazolo[4,5-b]pyridine"),
    # 1H-[1,2,3]triazolo[4,5-e]pyrazine (CH at 5)
    ("c1cnc2[nH]nnc2n1",        "1H-[1,2,3]triazolo[4,5-e]pyrazine"),
    ("Cc1cnc2[nH]nnc2n1",       "5-methyl-1H-[1,2,3]triazolo[4,5-e]pyrazine"),
    # 2H-[1,2,3]triazolo[4,5-e]pyrazine (CH at 5)
    ("c1cnc2n[nH]nc2n1",        "2H-[1,2,3]triazolo[4,5-e]pyrazine"),
    ("Cc1cnc2n[nH]nc2n1",       "5-methyl-2H-[1,2,3]triazolo[4,5-e]pyrazine"),
    # 1H-[1,2,3]triazolo[4,5-d]pyrimidine (CH at 5,7)
    ("c1ncc2[nH]nnc2n1",        "1H-[1,2,3]triazolo[4,5-d]pyrimidine"),
    ("Cc1ncc2[nH]nnc2n1",       "5-methyl-1H-[1,2,3]triazolo[4,5-d]pyrimidine"),
    ("Cc1ncnc2nn[nH]c12",       "7-methyl-1H-[1,2,3]triazolo[4,5-d]pyrimidine"),
    # 2H-[1,2,3]triazolo[4,5-d]pyrimidine (CH at 5,7)
    ("c1ncc2n[nH]nc2n1",        "2H-[1,2,3]triazolo[4,5-d]pyrimidine"),
    ("Cc1ncc2n[nH]nc2n1",       "5-methyl-2H-[1,2,3]triazolo[4,5-d]pyrimidine"),
    ("Cc1ncnc2n[nH]nc12",       "7-methyl-2H-[1,2,3]triazolo[4,5-d]pyrimidine"),
    # 1H-[1,2,3]triazolo[5,4-d]pyrimidine (CH at 5,7)
    ("c1ncc2nn[nH]c2n1",        "1H-[1,2,3]triazolo[5,4-d]pyrimidine"),
    ("Cc1ncc2nn[nH]c2n1",       "5-methyl-1H-[1,2,3]triazolo[5,4-d]pyrimidine"),
    ("Cc1ncnc2[nH]nnc12",       "7-methyl-1H-[1,2,3]triazolo[5,4-d]pyrimidine"),
    # 1H-[1,2,3]triazolo[4,5-e][1,2,4]triazine (CH at 6)
    ("c1nnc2[nH]nnc2n1",        "1H-[1,2,3]triazolo[4,5-e][1,2,4]triazine"),
    ("Cc1nnc2[nH]nnc2n1",       "6-methyl-1H-[1,2,3]triazolo[4,5-e][1,2,4]triazine"),
    # 2H-[1,2,3]triazolo[4,5-e][1,2,4]triazine (CH at 6)
    ("c1nnc2n[nH]nc2n1",        "2H-[1,2,3]triazolo[4,5-e][1,2,4]triazine"),
    ("Cc1nnc2n[nH]nc2n1",       "6-methyl-2H-[1,2,3]triazolo[4,5-e][1,2,4]triazine"),
    # 1H-[1,2,3]triazolo[5,4-e][1,2,4]triazine (CH at 6)
    ("c1nnc2nn[nH]c2n1",        "1H-[1,2,3]triazolo[5,4-e][1,2,4]triazine"),
    ("Cc1nnc2nn[nH]c2n1",       "6-methyl-1H-[1,2,3]triazolo[5,4-e][1,2,4]triazine"),
    # 1H-[1,2,3]triazolo[4,5-d][1,2,3]triazine (CH at 7)
    ("c1nnnc2nn[nH]c12",        "1H-[1,2,3]triazolo[4,5-d][1,2,3]triazine"),
    ("Cc1nnnc2nn[nH]c12",       "7-methyl-1H-[1,2,3]triazolo[4,5-d][1,2,3]triazine"),
    # 2H-[1,2,3]triazolo[4,5-d][1,2,3]triazine (CH at 7)
    ("c1nnnc2n[nH]nc12",        "2H-[1,2,3]triazolo[4,5-d][1,2,3]triazine"),
    ("Cc1nnnc2n[nH]nc12",       "7-methyl-2H-[1,2,3]triazolo[4,5-d][1,2,3]triazine"),
    # 1H-[1,2,3]triazolo[5,4-d][1,2,3]triazine (CH at 7)
    ("c1nnnc2[nH]nnc12",        "1H-[1,2,3]triazolo[5,4-d][1,2,3]triazine"),
    ("Cc1nnnc2[nH]nnc12",       "7-methyl-1H-[1,2,3]triazolo[5,4-d][1,2,3]triazine"),
    # 1H-pyrazolo[3,4-d][1,2,3]triazine (CH at 4,5)
    ("c1nnnc2n[nH]cc12",        "1H-pyrazolo[3,4-d][1,2,3]triazine"),
    ("Cc1nnnc2n[nH]cc12",       "4-methyl-1H-pyrazolo[3,4-d][1,2,3]triazine"),
    ("Cc1[nH]nc2nnncc12",       "5-methyl-1H-pyrazolo[3,4-d][1,2,3]triazine"),
    # 1H-pyrazolo[3,4-e][1,2,4]triazine (CH at 3,7)
    ("c1n[nH]c2cnnc-2n1",       "1H-pyrazolo[3,4-e][1,2,4]triazine"),
    ("Cc1n[nH]c2cnnc-2n1",      "3-methyl-1H-pyrazolo[3,4-e][1,2,4]triazine"),
    ("Cc1nnc2ncn[nH]c1-2",      "7-methyl-1H-pyrazolo[3,4-e][1,2,4]triazine"),
    # 1H-pyrazolo[4,3-d][1,2,3]triazine (CH at 4,7)
    ("c1nn[nH]c2cnnc1-2",       "1H-pyrazolo[4,3-d][1,2,3]triazine"),
    ("Cc1nn[nH]c2cnnc1-2",      "4-methyl-1H-pyrazolo[4,3-d][1,2,3]triazine"),
    ("Cc1nnc2cnn[nH]c1-2",      "7-methyl-1H-pyrazolo[4,3-d][1,2,3]triazine"),
    # 1H-pyrazolo[4,3-e][1,2,4]triazine (CH at 3,5)
    ("c1nnc2n[nH]cc2n1",        "1H-pyrazolo[4,3-e][1,2,4]triazine"),
    ("Cc1[nH]nc2nncnc12",       "3-methyl-1H-pyrazolo[4,3-e][1,2,4]triazine"),
    ("Cc1nnc2n[nH]cc2n1",       "5-methyl-1H-pyrazolo[4,3-e][1,2,4]triazine"),
    # 1H-pyrazolo[4,5-e][1,2,4]triazine (CH at 3,5)
    ("c1nnc2[nH]ncc2n1",        "1H-pyrazolo[4,5-e][1,2,4]triazine"),
    ("Cc1n[nH]c2nncnc12",       "3-methyl-1H-pyrazolo[4,5-e][1,2,4]triazine"),
    ("Cc1nnc2[nH]ncc2n1",       "5-methyl-1H-pyrazolo[4,5-e][1,2,4]triazine"),
    # 1H-pyrazolo[5,4-d][1,2,3]triazine (CH at 4,5)
    ("c1nnnc2[nH]ncc12",        "1H-pyrazolo[5,4-d][1,2,3]triazine"),
    ("Cc1nnnc2[nH]ncc12",       "4-methyl-1H-pyrazolo[5,4-d][1,2,3]triazine"),
    ("Cc1n[nH]c2nnncc12",       "5-methyl-1H-pyrazolo[5,4-d][1,2,3]triazine"),
    # 1,2,4-triazine (CH at 3,6)
    ("c1cnncn1",                "1,2,4-triazine"),
    ("Cc1nccnn1",               "3-methyl-1,2,4-triazine"),
    ("Cc1cncnn1",               "6-methyl-1,2,4-triazine"),
    # 1,2-thiazole (isothiazole PIN)
    ("c1cnsc1",                 "1,2-thiazole"),
    ("Cc1ccsn1",                "3-methyl-1,2-thiazole"),
    # 1,3-oxazole (oxazole PIN)
    ("c1cocn1",                 "1,3-oxazole"),
    ("Cc1ncco1",                "2-methyl-1,3-oxazole"),
    # 1,3-thiazole (thiazole PIN)
    ("c1cscn1",                 "1,3-thiazole"),
    ("Cc1nccs1",                "2-methyl-1,3-thiazole"),
])
def test_phase729(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
