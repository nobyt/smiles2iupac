"""Phase 726: missing methyl derivatives of partially saturated benzo-fused
heterocycles: 1,3-dihydro-2-benzothiophene, 2,3-dihydrobenzothiophene,
2,1-benzisothiazole, 3,4-dihydroisoquinoline, 3,4-dihydroquinoline,
3,4-dihydro-2H-1,4-benzothiazine, thiochroman,
2,1,3-benzothiadiazole, 2,1,3-benzoxadiazole, 2,1,3-benzoselenadiazole.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,3-dihydro-2-benzothiophene (CH at 1,4,5; C2v: 1=3, 4=7, 5=6)
    ("c1ccc2c(c1)CSC2",       "1,3-dihydro-2-benzothiophene"),
    ("CC1SCc2ccccc21",        "1-methyl-1,3-dihydro-2-benzothiophene"),
    ("Cc1cccc2c1CSC2",        "4-methyl-1,3-dihydro-2-benzothiophene"),
    ("Cc1ccc2c(c1)CSC2",      "5-methyl-1,3-dihydro-2-benzothiophene"),
    # 2,3-dihydrobenzothiophene (CH at 2,3,4,5,6,7; no symmetry)
    ("c1ccc2c(c1)CCS2",       "2,3-dihydrobenzothiophene"),
    ("CC1Cc2ccccc2S1",        "2-methyl-2,3-dihydrobenzothiophene"),
    ("CC1CSc2ccccc21",        "3-methyl-2,3-dihydrobenzothiophene"),
    # 2,1-benzisothiazole (CH at 3,4,5,6,7)
    ("c1ccc2nscc2c1",         "2,1-benzisothiazole"),
    ("Cc1snc2ccccc12",        "3-methyl-2,1-benzisothiazole"),
    ("Cc1cccc2nscc12",        "4-methyl-2,1-benzisothiazole"),
    ("Cc1ccc2csnc2c1",        "6-methyl-2,1-benzisothiazole"),
    ("Cc1cccc2csnc12",        "7-methyl-2,1-benzisothiazole"),
    # 3,4-dihydroisoquinoline (CH at 1,3,4,5,6,7,8)
    ("C1=NCCc2ccccc21",       "3,4-dihydroisoquinoline"),
    ("Cc1ccc2c(c1)CCN=C2",    "6-methyl-3,4-dihydroisoquinoline"),
    # 3,4-dihydroquinoline (CH at 2,3,4,5,6,7,8)
    ("C1=Nc2ccccc2CC1",       "3,4-dihydroquinoline"),
    ("Cc1ccc2c(c1)CCC=N2",    "6-methyl-3,4-dihydroquinoline"),
    ("Cc1ccc2c(c1)N=CCC2",    "7-methyl-3,4-dihydroquinoline"),
    # 3,4-dihydro-2H-1,4-benzothiazine (CH at 2,3,5,6,7,8)
    ("c1ccc2c(c1)NCCS2",      "3,4-dihydro-2H-1,4-benzothiazine"),
    ("CC1CNc2ccccc2S1",       "2-methyl-3,4-dihydro-2H-1,4-benzothiazine"),
    ("CC1CSc2ccccc2N1",       "3-methyl-3,4-dihydro-2H-1,4-benzothiazine"),
    ("Cc1cccc2c1NCCS2",       "5-methyl-3,4-dihydro-2H-1,4-benzothiazine"),
    ("Cc1ccc2c(c1)NCCS2",     "6-methyl-3,4-dihydro-2H-1,4-benzothiazine"),
    ("Cc1ccc2c(c1)SCCN2",     "7-methyl-3,4-dihydro-2H-1,4-benzothiazine"),
    ("Cc1cccc2c1SCCN2",       "8-methyl-3,4-dihydro-2H-1,4-benzothiazine"),
    # thiochroman (CH at 2,3,4,5,6,7,8)
    ("c1ccc2c(c1)CCCS2",      "thiochroman"),
    ("Cc1ccc2c(c1)SCCC2",     "7-methylthiochroman"),
    ("Cc1cccc2c1SCCC2",       "8-methylthiochroman"),
    # 2,1,3-benzothiadiazole (CH at 4,5; C2v: 4=7, 5=6)
    ("c1ccc2nsnc2c1",         "2,1,3-benzothiadiazole"),
    ("Cc1cccc2nsnc12",        "4-methyl-2,1,3-benzothiadiazole"),
    # 2,1,3-benzoxadiazole (CH at 4,5; C2v: 4=7, 5=6)
    ("c1ccc2nonc2c1",         "2,1,3-benzoxadiazole"),
    ("Cc1cccc2nonc12",        "4-methyl-2,1,3-benzoxadiazole"),
    # 2,1,3-benzoselenadiazole (CH at 4,5; C2v: 4=7, 5=6)
    ("c1ccc2n[se]nc2c1",      "2,1,3-benzoselenadiazole"),
    ("Cc1cccc2n[se]nc12",     "4-methyl-2,1,3-benzoselenadiazole"),
])
def test_phase726(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
