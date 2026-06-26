"""Phase 725: methyl derivatives of partially saturated benzo-fused heterocycles:
2,3-dihydro-1H-benzimidazole, 2,3-dihydro-1H-indazole,
2,3-dihydrobenzo[d]oxazole, 2,3-dihydrobenzo[d]thiazole,
2,3-dihydrobenzo[d]isothiazole, 2,3-dihydrobenzo[d]isoxazole,
4,5,6,7-tetrahydrobenzofuran, 4,5,6,7-tetrahydrobenzothiophene,
4,5,6,7-tetrahydrobenzo[d]oxazole, 4,5,6,7-tetrahydrobenzo[d]thiazole,
4,5,6,7-tetrahydrobenzo[d]isothiazole, 4,5,6,7-tetrahydrobenzo[d]isoxazole,
4,5,6,7-tetrahydro-1H-benzimidazole, 4,5,6,7-tetrahydro-1H-indazole,
4,5,6,7-tetrahydro-1H-indole,
4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole,
4,5,6,7-tetrahydro-2H-benzo[d][1,2,3]triazole.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2,3-dihydro-1H-benzimidazole (CH at 2,4,5; C2v: 4≡7, 5≡6)
    ("c1ccc2c(c1)NCN2",                                   "2,3-dihydro-1H-benzimidazole"),
    ("CC1Nc2ccccc2N1",                                    "2-methyl-2,3-dihydro-1H-benzimidazole"),
    ("Cc1cccc2c1NCN2",                                    "4-methyl-2,3-dihydro-1H-benzimidazole"),
    ("Cc1ccc2c(c1)NCN2",                                  "5-methyl-2,3-dihydro-1H-benzimidazole"),
    # 2,3-dihydro-1H-indazole (CH at 3,4,5,6,7)
    ("c1ccc2c(c1)CNN2",                                   "2,3-dihydro-1H-indazole"),
    ("CC1NNc2ccccc21",                                    "3-methyl-2,3-dihydro-1H-indazole"),
    ("Cc1cccc2c1CNN2",                                    "4-methyl-2,3-dihydro-1H-indazole"),
    ("Cc1ccc2c(c1)CNN2",                                  "5-methyl-2,3-dihydro-1H-indazole"),
    ("Cc1ccc2c(c1)NNC2",                                  "6-methyl-2,3-dihydro-1H-indazole"),
    ("Cc1cccc2c1NNC2",                                    "7-methyl-2,3-dihydro-1H-indazole"),
    # 2,3-dihydrobenzo[d]oxazole (CH at 2,4,5,6,7)
    ("c1ccc2c(c1)NCO2",                                   "2,3-dihydrobenzo[d]oxazole"),
    ("CC1Nc2ccccc2O1",                                    "2-methyl-2,3-dihydrobenzo[d]oxazole"),
    ("Cc1cccc2c1NCO2",                                    "4-methyl-2,3-dihydrobenzo[d]oxazole"),
    ("Cc1ccc2c(c1)NCO2",                                  "5-methyl-2,3-dihydrobenzo[d]oxazole"),
    ("Cc1ccc2c(c1)OCN2",                                  "6-methyl-2,3-dihydrobenzo[d]oxazole"),
    ("Cc1cccc2c1OCN2",                                    "7-methyl-2,3-dihydrobenzo[d]oxazole"),
    # 2,3-dihydrobenzo[d]thiazole (CH at 2,4,5,6,7)
    ("c1ccc2c(c1)NCS2",                                   "2,3-dihydrobenzo[d]thiazole"),
    ("CC1Nc2ccccc2S1",                                    "2-methyl-2,3-dihydrobenzo[d]thiazole"),
    ("Cc1cccc2c1NCS2",                                    "4-methyl-2,3-dihydrobenzo[d]thiazole"),
    ("Cc1ccc2c(c1)NCS2",                                  "5-methyl-2,3-dihydrobenzo[d]thiazole"),
    ("Cc1ccc2c(c1)SCN2",                                  "6-methyl-2,3-dihydrobenzo[d]thiazole"),
    ("Cc1cccc2c1SCN2",                                    "7-methyl-2,3-dihydrobenzo[d]thiazole"),
    # 2,3-dihydrobenzo[d]isothiazole (CH at 3,4,5,6,7)
    ("c1ccc2c(c1)CNS2",                                   "2,3-dihydrobenzo[d]isothiazole"),
    ("CC1NSc2ccccc21",                                    "3-methyl-2,3-dihydrobenzo[d]isothiazole"),
    ("Cc1cccc2c1CNS2",                                    "4-methyl-2,3-dihydrobenzo[d]isothiazole"),
    ("Cc1ccc2c(c1)CNS2",                                  "5-methyl-2,3-dihydrobenzo[d]isothiazole"),
    ("Cc1ccc2c(c1)SNC2",                                  "6-methyl-2,3-dihydrobenzo[d]isothiazole"),
    ("Cc1cccc2c1SNC2",                                    "7-methyl-2,3-dihydrobenzo[d]isothiazole"),
    # 2,3-dihydrobenzo[d]isoxazole (CH at 3,4,5,6,7)
    ("c1ccc2c(c1)CNO2",                                   "2,3-dihydrobenzo[d]isoxazole"),
    ("CC1NOc2ccccc21",                                    "3-methyl-2,3-dihydrobenzo[d]isoxazole"),
    ("Cc1cccc2c1CNO2",                                    "4-methyl-2,3-dihydrobenzo[d]isoxazole"),
    ("Cc1ccc2c(c1)CNO2",                                  "5-methyl-2,3-dihydrobenzo[d]isoxazole"),
    ("Cc1ccc2c(c1)ONC2",                                  "6-methyl-2,3-dihydrobenzo[d]isoxazole"),
    ("Cc1cccc2c1ONC2",                                    "7-methyl-2,3-dihydrobenzo[d]isoxazole"),
    # 4,5,6,7-tetrahydrobenzofuran (CH at 2,3,4,5,6,7)
    ("c1cc2c(o1)CCCC2",                                   "4,5,6,7-tetrahydrobenzofuran"),
    ("Cc1cc2c(o1)CCCC2",                                  "2-methyl-4,5,6,7-tetrahydrobenzofuran"),
    ("Cc1coc2c1CCCC2",                                    "3-methyl-4,5,6,7-tetrahydrobenzofuran"),
    ("CC1CCCc2occc21",                                    "4-methyl-4,5,6,7-tetrahydrobenzofuran"),
    ("CC1CCc2occc2C1",                                    "5-methyl-4,5,6,7-tetrahydrobenzofuran"),
    ("CC1CCc2ccoc2C1",                                    "6-methyl-4,5,6,7-tetrahydrobenzofuran"),
    ("CC1CCCc2ccoc21",                                    "7-methyl-4,5,6,7-tetrahydrobenzofuran"),
    # 4,5,6,7-tetrahydrobenzothiophene (CH at 2,3,4,5,6,7)
    ("c1cc2c(s1)CCCC2",                                   "4,5,6,7-tetrahydrobenzothiophene"),
    ("Cc1cc2c(s1)CCCC2",                                  "2-methyl-4,5,6,7-tetrahydrobenzothiophene"),
    ("Cc1csc2c1CCCC2",                                    "3-methyl-4,5,6,7-tetrahydrobenzothiophene"),
    ("CC1CCCc2sccc21",                                    "4-methyl-4,5,6,7-tetrahydrobenzothiophene"),
    ("CC1CCc2sccc2C1",                                    "5-methyl-4,5,6,7-tetrahydrobenzothiophene"),
    ("CC1CCc2ccsc2C1",                                    "6-methyl-4,5,6,7-tetrahydrobenzothiophene"),
    ("CC1CCCc2ccsc21",                                    "7-methyl-4,5,6,7-tetrahydrobenzothiophene"),
    # 4,5,6,7-tetrahydrobenzo[d]oxazole (CH at 2,4,5,6,7)
    ("c1nc2c(o1)CCCC2",                                   "4,5,6,7-tetrahydrobenzo[d]oxazole"),
    ("Cc1nc2c(o1)CCCC2",                                  "2-methyl-4,5,6,7-tetrahydrobenzo[d]oxazole"),
    ("CC1CCCc2ocnc21",                                    "4-methyl-4,5,6,7-tetrahydrobenzo[d]oxazole"),
    ("CC1CCc2ocnc2C1",                                    "5-methyl-4,5,6,7-tetrahydrobenzo[d]oxazole"),
    ("CC1CCc2ncoc2C1",                                    "6-methyl-4,5,6,7-tetrahydrobenzo[d]oxazole"),
    ("CC1CCCc2ncoc21",                                    "7-methyl-4,5,6,7-tetrahydrobenzo[d]oxazole"),
    # 4,5,6,7-tetrahydrobenzo[d]thiazole (CH at 2,4,5,6,7)
    ("c1nc2c(s1)CCCC2",                                   "4,5,6,7-tetrahydrobenzo[d]thiazole"),
    ("Cc1nc2c(s1)CCCC2",                                  "2-methyl-4,5,6,7-tetrahydrobenzo[d]thiazole"),
    ("CC1CCCc2scnc21",                                    "4-methyl-4,5,6,7-tetrahydrobenzo[d]thiazole"),
    ("CC1CCc2scnc2C1",                                    "5-methyl-4,5,6,7-tetrahydrobenzo[d]thiazole"),
    ("CC1CCc2ncsc2C1",                                    "6-methyl-4,5,6,7-tetrahydrobenzo[d]thiazole"),
    ("CC1CCCc2ncsc21",                                    "7-methyl-4,5,6,7-tetrahydrobenzo[d]thiazole"),
    # 4,5,6,7-tetrahydrobenzo[d]isothiazole (CH at 3,4,5,6,7)
    ("c1nsc2c1CCCC2",                                     "4,5,6,7-tetrahydrobenzo[d]isothiazole"),
    ("Cc1nsc2c1CCCC2",                                    "3-methyl-4,5,6,7-tetrahydrobenzo[d]isothiazole"),
    ("CC1CCCc2sncc21",                                    "4-methyl-4,5,6,7-tetrahydrobenzo[d]isothiazole"),
    ("CC1CCc2sncc2C1",                                    "5-methyl-4,5,6,7-tetrahydrobenzo[d]isothiazole"),
    ("CC1CCc2cnsc2C1",                                    "6-methyl-4,5,6,7-tetrahydrobenzo[d]isothiazole"),
    ("CC1CCCc2cnsc21",                                    "7-methyl-4,5,6,7-tetrahydrobenzo[d]isothiazole"),
    # 4,5,6,7-tetrahydrobenzo[d]isoxazole (CH at 3,4,5,6,7)
    ("c1noc2c1CCCC2",                                     "4,5,6,7-tetrahydrobenzo[d]isoxazole"),
    ("Cc1noc2c1CCCC2",                                    "3-methyl-4,5,6,7-tetrahydrobenzo[d]isoxazole"),
    ("CC1CCCc2oncc21",                                    "4-methyl-4,5,6,7-tetrahydrobenzo[d]isoxazole"),
    ("CC1CCc2oncc2C1",                                    "5-methyl-4,5,6,7-tetrahydrobenzo[d]isoxazole"),
    ("CC1CCc2cnoc2C1",                                    "6-methyl-4,5,6,7-tetrahydrobenzo[d]isoxazole"),
    ("CC1CCCc2cnoc21",                                    "7-methyl-4,5,6,7-tetrahydrobenzo[d]isoxazole"),
    # 4,5,6,7-tetrahydro-1H-benzimidazole (CH at 2,4,5; C2v: 4≡7, 5≡6)
    ("c1nc2c([nH]1)CCCC2",                                "4,5,6,7-tetrahydro-1H-benzimidazole"),
    ("Cc1nc2c([nH]1)CCCC2",                               "2-methyl-4,5,6,7-tetrahydro-1H-benzimidazole"),
    ("CC1CCCc2nc[nH]c21",                                 "4-methyl-4,5,6,7-tetrahydro-1H-benzimidazole"),
    ("CC1CCc2nc[nH]c2C1",                                 "5-methyl-4,5,6,7-tetrahydro-1H-benzimidazole"),
    # 4,5,6,7-tetrahydro-1H-indazole (CH at 3,4,5,6,7)
    ("c1n[nH]c2c1CCCC2",                                  "4,5,6,7-tetrahydro-1H-indazole"),
    ("Cc1n[nH]c2c1CCCC2",                                 "3-methyl-4,5,6,7-tetrahydro-1H-indazole"),
    ("CC1CCCc2[nH]ncc21",                                 "4-methyl-4,5,6,7-tetrahydro-1H-indazole"),
    ("CC1CCc2[nH]ncc2C1",                                 "5-methyl-4,5,6,7-tetrahydro-1H-indazole"),
    ("CC1CCc2cn[nH]c2C1",                                 "6-methyl-4,5,6,7-tetrahydro-1H-indazole"),
    ("CC1CCCc2cn[nH]c21",                                 "7-methyl-4,5,6,7-tetrahydro-1H-indazole"),
    # 4,5,6,7-tetrahydro-1H-indole (CH at 2,3,4,5,6,7)
    ("c1cc2c([nH]1)CCCC2",                                "4,5,6,7-tetrahydro-1H-indole"),
    ("Cc1cc2c([nH]1)CCCC2",                               "2-methyl-4,5,6,7-tetrahydro-1H-indole"),
    ("Cc1c[nH]c2c1CCCC2",                                 "3-methyl-4,5,6,7-tetrahydro-1H-indole"),
    ("CC1CCCc2[nH]ccc21",                                 "4-methyl-4,5,6,7-tetrahydro-1H-indole"),
    ("CC1CCc2[nH]ccc2C1",                                 "5-methyl-4,5,6,7-tetrahydro-1H-indole"),
    ("CC1CCc2cc[nH]c2C1",                                 "6-methyl-4,5,6,7-tetrahydro-1H-indole"),
    ("CC1CCCc2cc[nH]c21",                                 "7-methyl-4,5,6,7-tetrahydro-1H-indole"),
    # 4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole (CH at 4,5,6,7)
    ("C1CCc2[nH]nnc2C1",                                  "4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole"),
    ("CC1CCCc2[nH]nnc21",                                 "4-methyl-4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole"),
    ("CC1CCc2[nH]nnc2C1",                                 "5-methyl-4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole"),
    ("CC1CCc2nn[nH]c2C1",                                 "6-methyl-4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole"),
    ("CC1CCCc2nn[nH]c21",                                 "7-methyl-4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole"),
    # 4,5,6,7-tetrahydro-2H-benzo[d][1,2,3]triazole (CH at 4,5; C2v: 4≡7, 5≡6)
    ("C1CCc2n[nH]nc2C1",                                  "4,5,6,7-tetrahydro-2H-benzo[d][1,2,3]triazole"),
    ("CC1CCCc2n[nH]nc21",                                 "4-methyl-4,5,6,7-tetrahydro-2H-benzo[d][1,2,3]triazole"),
    ("CC1CCc2n[nH]nc2C1",                                 "5-methyl-4,5,6,7-tetrahydro-2H-benzo[d][1,2,3]triazole"),
])
def test_phase725(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
