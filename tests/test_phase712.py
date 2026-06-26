"""Phase 712: methyl derivatives of simple benzene-fused aromatic heterocycles
(benzofuran, benzo[b]thiophene, benzo[b]selenophene, 1,3-benzodioxole,
1,3-benzoxazole, 1,3-benzothiazole, 1,2-benzisoxazole, 1,2-benzisothiazole,
1H-benzimidazole, 1H-indole, 1H-indazole, 2H-indazole,
1H-benzotriazole, 1,2,3-benzothiadiazole, 1,2,3-benzoselenadiazole).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzofuran (CH at 2,3,4,5,6,7)
    ("c1ccc2occc2c1",       "benzofuran"),
    ("Cc1cc2ccccc2o1",      "2-methylbenzofuran"),
    ("Cc1coc2ccccc12",      "3-methylbenzofuran"),
    ("Cc1cccc2occc12",      "4-methylbenzofuran"),
    ("Cc1ccc2occc2c1",      "5-methylbenzofuran"),
    ("Cc1ccc2ccoc2c1",      "6-methylbenzofuran"),
    ("Cc1cccc2ccoc12",      "7-methylbenzofuran"),
    # benzo[b]thiophene (CH at 2,3,4,5,6,7)
    ("c1ccc2sccc2c1",       "benzo[b]thiophene"),
    ("Cc1cc2ccccc2s1",      "2-methylbenzo[b]thiophene"),
    ("Cc1csc2ccccc12",      "3-methylbenzo[b]thiophene"),
    ("Cc1cccc2sccc12",      "4-methylbenzo[b]thiophene"),
    ("Cc1ccc2sccc2c1",      "5-methylbenzo[b]thiophene"),
    ("Cc1ccc2ccsc2c1",      "6-methylbenzo[b]thiophene"),
    ("Cc1cccc2ccsc12",      "7-methylbenzo[b]thiophene"),
    # benzo[b]selenophene (CH at 2,3,4,5,6,7)
    ("c1ccc2[se]ccc2c1",    "benzo[b]selenophene"),
    ("Cc1cc2ccccc2[se]1",   "2-methylbenzo[b]selenophene"),
    ("Cc1c[se]c2ccccc12",   "3-methylbenzo[b]selenophene"),
    ("Cc1cccc2[se]ccc12",   "4-methylbenzo[b]selenophene"),
    ("Cc1ccc2[se]ccc2c1",   "5-methylbenzo[b]selenophene"),
    ("Cc1ccc2cc[se]c2c1",   "6-methylbenzo[b]selenophene"),
    ("Cc1cccc2cc[se]c12",   "7-methylbenzo[b]selenophene"),
    # 1,3-benzodioxole (CH at 4,5)
    ("c1ccc2c(c1)OCO2",     "1,3-benzodioxole"),
    ("Cc1cccc2c1OCO2",      "4-methyl-1,3-benzodioxole"),
    ("Cc1ccc2c(c1)OCO2",    "5-methyl-1,3-benzodioxole"),
    # 1,3-benzoxazole (CH at 2,4,5,6,7)
    ("c1ccc2ocnc2c1",       "1,3-benzoxazole"),
    ("Cc1nc2ccccc2o1",      "2-methyl-1,3-benzoxazole"),
    ("Cc1cccc2ocnc12",      "4-methyl-1,3-benzoxazole"),
    ("Cc1ccc2ocnc2c1",      "5-methyl-1,3-benzoxazole"),
    ("Cc1ccc2ncoc2c1",      "6-methyl-1,3-benzoxazole"),
    ("Cc1cccc2ncoc12",      "7-methyl-1,3-benzoxazole"),
    # 1,3-benzothiazole (CH at 2,4,5,6,7)
    ("c1ccc2scnc2c1",       "1,3-benzothiazole"),
    ("Cc1nc2ccccc2s1",      "2-methyl-1,3-benzothiazole"),
    ("Cc1cccc2scnc12",      "4-methyl-1,3-benzothiazole"),
    ("Cc1ccc2scnc2c1",      "5-methyl-1,3-benzothiazole"),
    ("Cc1ccc2ncsc2c1",      "6-methyl-1,3-benzothiazole"),
    ("Cc1cccc2ncsc12",      "7-methyl-1,3-benzothiazole"),
    # 1,2-benzisoxazole (CH at 3,4,5,6,7)
    ("c1ccc2oncc2c1",       "1,2-benzisoxazole"),
    ("Cc1noc2ccccc12",      "3-methyl-1,2-benzisoxazole"),
    ("Cc1cccc2oncc12",      "4-methyl-1,2-benzisoxazole"),
    ("Cc1ccc2oncc2c1",      "5-methyl-1,2-benzisoxazole"),
    ("Cc1ccc2cnoc2c1",      "6-methyl-1,2-benzisoxazole"),
    ("Cc1cccc2cnoc12",      "7-methyl-1,2-benzisoxazole"),
    # 1,2-benzisothiazole (CH at 3,4,5,6,7)
    ("c1ccc2sncc2c1",       "1,2-benzisothiazole"),
    ("Cc1nsc2ccccc12",      "3-methyl-1,2-benzisothiazole"),
    ("Cc1cccc2sncc12",      "4-methyl-1,2-benzisothiazole"),
    ("Cc1ccc2sncc2c1",      "5-methyl-1,2-benzisothiazole"),
    ("Cc1ccc2cnsc2c1",      "6-methyl-1,2-benzisothiazole"),
    ("Cc1cccc2cnsc12",      "7-methyl-1,2-benzisothiazole"),
    # 1H-benzimidazole (CH at 2,4,5,6,7)
    ("c1ccc2[nH]cnc2c1",    "1H-benzimidazole"),
    ("Cc1nc2ccccc2[nH]1",   "2-methyl-1H-benzimidazole"),
    ("Cc1cccc2[nH]cnc12",   "4-methyl-1H-benzimidazole"),
    ("Cc1ccc2[nH]cnc2c1",   "5-methyl-1H-benzimidazole"),
    ("Cc1ccc2nc[nH]c2c1",   "6-methyl-1H-benzimidazole"),
    ("Cc1cccc2nc[nH]c12",   "7-methyl-1H-benzimidazole"),
    # 1H-indole (CH at 2,3,4,5,6,7)
    ("c1ccc2[nH]ccc2c1",    "1H-indole"),
    ("Cc1cc2ccccc2[nH]1",   "2-methyl-1H-indole"),
    ("Cc1c[nH]c2ccccc12",   "3-methyl-1H-indole"),
    ("Cc1cccc2[nH]ccc12",   "4-methyl-1H-indole"),
    ("Cc1ccc2[nH]ccc2c1",   "5-methyl-1H-indole"),
    ("Cc1ccc2cc[nH]c2c1",   "6-methyl-1H-indole"),
    ("Cc1cccc2cc[nH]c12",   "7-methyl-1H-indole"),
    # 1H-indazole (CH at 3,4,5,6,7)
    ("c1ccc2[nH]ncc2c1",    "1H-indazole"),
    ("Cc1n[nH]c2ccccc12",   "3-methyl-1H-indazole"),
    ("Cc1cccc2[nH]ncc12",   "4-methyl-1H-indazole"),
    ("Cc1ccc2[nH]ncc2c1",   "5-methyl-1H-indazole"),
    ("Cc1ccc2cn[nH]c2c1",   "6-methyl-1H-indazole"),
    ("Cc1cccc2cn[nH]c12",   "7-methyl-1H-indazole"),
    # 2H-indazole (CH at 3,4,5,6,7)
    ("c1ccc2n[nH]cc2c1",    "2H-indazole"),
    ("Cc1[nH]nc2ccccc12",   "3-methyl-2H-indazole"),
    ("Cc1cccc2n[nH]cc12",   "4-methyl-2H-indazole"),
    ("Cc1ccc2n[nH]cc2c1",   "5-methyl-2H-indazole"),
    ("Cc1ccc2c[nH]nc2c1",   "6-methyl-2H-indazole"),
    ("Cc1cccc2c[nH]nc12",   "7-methyl-2H-indazole"),
    # 1H-benzotriazole (CH at 4,5,6,7)
    ("c1ccc2[nH]nnc2c1",    "1H-benzotriazole"),
    ("Cc1cccc2[nH]nnc12",   "4-methyl-1H-benzotriazole"),
    ("Cc1ccc2[nH]nnc2c1",   "5-methyl-1H-benzotriazole"),
    ("Cc1ccc2nn[nH]c2c1",   "6-methyl-1H-benzotriazole"),
    ("Cc1cccc2nn[nH]c12",   "7-methyl-1H-benzotriazole"),
    # 1,2,3-benzothiadiazole (CH at 4,5,6,7)
    ("c1ccc2snnc2c1",       "1,2,3-benzothiadiazole"),
    ("Cc1cccc2snnc12",      "4-methyl-1,2,3-benzothiadiazole"),
    ("Cc1ccc2snnc2c1",      "5-methyl-1,2,3-benzothiadiazole"),
    ("Cc1ccc2nnsc2c1",      "6-methyl-1,2,3-benzothiadiazole"),
    ("Cc1cccc2nnsc12",      "7-methyl-1,2,3-benzothiadiazole"),
    # 1,2,3-benzoselenadiazole (CH at 4,5,6,7)
    ("c1ccc2[se]nnc2c1",    "1,2,3-benzoselenadiazole"),
    ("Cc1cccc2[se]nnc12",   "4-methyl-1,2,3-benzoselenadiazole"),
    ("Cc1ccc2[se]nnc2c1",   "5-methyl-1,2,3-benzoselenadiazole"),
    ("Cc1ccc2nn[se]c2c1",   "6-methyl-1,2,3-benzoselenadiazole"),
    ("Cc1cccc2nn[se]c12",   "7-methyl-1,2,3-benzoselenadiazole"),
])
def test_phase712(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
