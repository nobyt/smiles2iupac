"""Phase 622: methyl-substituted benzo-fused heterocycles (IUPAC 2013).

All carbon positions on quinoline, isoquinoline, 1H-indole, coumarin,
benzofuran, benzo[b]thiophene, 1H-benzimidazole, 1,3-benzothiazole,
and 1,3-benzoxazole.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # quinoline (N at 1)
    ("Cc1ccc2ncccc2c1",   "6-methylquinoline"),
    ("Cc1ccc2cccnc2c1",   "7-methylquinoline"),
    ("Cc1cccc2cccnc12",   "8-methylquinoline"),
    ("Cc1ccc2ccccc2n1",   "2-methylquinoline"),
    ("Cc1cnc2ccccc2c1",   "3-methylquinoline"),
    ("Cc1ccnc2ccccc12",   "4-methylquinoline"),
    ("Cc1cccc2ncccc12",   "5-methylquinoline"),
    # isoquinoline (N at 2)
    ("Cc1ccc2cnccc2c1",   "6-methylisoquinoline"),
    ("Cc1ccc2ccncc2c1",   "7-methylisoquinoline"),
    ("Cc1cccc2ccncc12",   "8-methylisoquinoline"),
    ("Cc1nccc2ccccc12",   "1-methylisoquinoline"),
    ("Cc1cc2ccccc2cn1",   "3-methylisoquinoline"),
    ("Cc1cncc2ccccc12",   "4-methylisoquinoline"),
    ("Cc1cccc2cnccc12",   "5-methylisoquinoline"),
    # 1H-indole (NH at 1, C at 2–7)
    ("Cc1ccc2[nH]ccc2c1", "5-methyl-1H-indole"),
    ("Cc1ccc2cc[nH]c2c1", "6-methyl-1H-indole"),
    ("Cc1cccc2cc[nH]c12", "7-methyl-1H-indole"),
    ("Cc1cc2ccccc2[nH]1", "2-methyl-1H-indole"),
    ("Cc1c[nH]c2ccccc12", "3-methyl-1H-indole"),
    ("Cc1cccc2[nH]ccc12", "4-methyl-1H-indole"),
    # coumarin (O at 1, C=O at 2)
    ("Cc1cc2ccccc2oc1=O",  "3-methylcoumarin"),
    ("Cc1cc(=O)oc2ccccc12", "4-methylcoumarin"),
    ("Cc1cccc2oc(=O)ccc12", "5-methylcoumarin"),
    ("Cc1ccc2oc(=O)ccc2c1", "6-methylcoumarin"),
    ("Cc1ccc2ccc(=O)oc2c1", "7-methylcoumarin"),
    ("Cc1cccc2ccc(=O)oc12", "8-methylcoumarin"),
    # benzofuran (O at 1)
    ("Cc1ccc2occc2c1",   "5-methylbenzofuran"),
    ("Cc1ccc2ccoc2c1",   "6-methylbenzofuran"),
    ("Cc1cccc2ccoc12",   "7-methylbenzofuran"),
    ("Cc1cc2ccccc2o1",   "2-methylbenzofuran"),
    ("Cc1coc2ccccc12",   "3-methylbenzofuran"),
    ("Cc1cccc2occc12",   "4-methylbenzofuran"),
    # benzo[b]thiophene (S at 1)
    ("Cc1ccc2sccc2c1",   "5-methylbenzo[b]thiophene"),
    ("Cc1ccc2ccsc2c1",   "6-methylbenzo[b]thiophene"),
    ("Cc1cccc2ccsc12",   "7-methylbenzo[b]thiophene"),
    ("Cc1cc2ccccc2s1",   "2-methylbenzo[b]thiophene"),
    ("Cc1csc2ccccc12",   "3-methylbenzo[b]thiophene"),
    ("Cc1cccc2sccc12",   "4-methylbenzo[b]thiophene"),
    # 1H-benzimidazole (NH at 1, N at 3)
    ("Cc1ccc2[nH]cnc2c1", "5-methyl-1H-benzimidazole"),
    ("Cc1ccc2nc[nH]c2c1", "6-methyl-1H-benzimidazole"),
    ("Cc1cccc2nc[nH]c12", "7-methyl-1H-benzimidazole"),
    ("Cc1nc2ccccc2[nH]1", "2-methyl-1H-benzimidazole"),
    ("Cc1cccc2[nH]cnc12", "4-methyl-1H-benzimidazole"),
    # 1,3-benzothiazole (S at 1, N at 3)
    ("Cc1ccc2scnc2c1",   "5-methyl-1,3-benzothiazole"),
    ("Cc1ccc2ncsc2c1",   "6-methyl-1,3-benzothiazole"),
    ("Cc1cccc2ncsc12",   "7-methyl-1,3-benzothiazole"),
    ("Cc1nc2ccccc2s1",   "2-methyl-1,3-benzothiazole"),
    ("Cc1cccc2scnc12",   "4-methyl-1,3-benzothiazole"),
    # 1,3-benzoxazole (O at 1, N at 3)
    ("Cc1ccc2ocnc2c1",   "5-methyl-1,3-benzoxazole"),
    ("Cc1ccc2ncoc2c1",   "6-methyl-1,3-benzoxazole"),
    ("Cc1cccc2ncoc12",   "7-methyl-1,3-benzoxazole"),
    ("Cc1nc2ccccc2o1",   "2-methyl-1,3-benzoxazole"),
    ("Cc1cccc2ocnc12",   "4-methyl-1,3-benzoxazole"),
])
def test_phase622_methyl_benzo_fused_heterocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
