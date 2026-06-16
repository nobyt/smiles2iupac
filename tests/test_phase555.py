"""Phase 555: Substituted benzo-fused thiadiazoles, oxadiazole, and selenadiazoles.
Only benzo-ring carbons (C4-C7) are substitutable in these fully heteroatom
5-membered rings fused to benzene.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2,3-benzothiadiazole
    ("c1ccc2snnc2c1",        "1,2,3-benzothiadiazole"),
    ("Cc1ccc2snnc2c1",       "5-methyl-1,2,3-benzothiadiazole"),
    ("Cc1cccc2snnc12",       "4-methyl-1,2,3-benzothiadiazole"),
    # 2,1,3-benzothiadiazole
    ("c1ccc2nsnc2c1",        "2,1,3-benzothiadiazole"),
    ("Cc1ccc2nsnc2c1",       "5-methyl-2,1,3-benzothiadiazole"),
    # 2,1,3-benzoxadiazole
    ("c1ccc2nonc2c1",        "2,1,3-benzoxadiazole"),
    ("Cc1ccc2nonc2c1",       "5-methyl-2,1,3-benzoxadiazole"),
    # 2,1,3-benzoselenadiazole
    ("c1ccc2n[se]nc2c1",     "2,1,3-benzoselenadiazole"),
    ("Cc1ccc2n[se]nc2c1",    "5-methyl-2,1,3-benzoselenadiazole"),
    # 1,2,3-benzoselenadiazole
    ("c1ccc2[se]nnc2c1",     "1,2,3-benzoselenadiazole"),
])
def test_phase555_benzo_diazoles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
