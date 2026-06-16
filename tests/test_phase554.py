"""Phase 554: Substituted 1,3-benzothiazole, 1,3-benzoxazole,
1,2-benzisothiazole, 1,2-benzisoxazole, and 2,1-benzisothiazole.
All share 9-atom bicyclic topology; substituent-capable carbons differ.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,3-benzothiazole: substitutable at C2, C4, C5, C6, C7
    ("c1ccc2scnc2c1",        "1,3-benzothiazole"),
    ("Cc1nc2ccccc2s1",       "2-methyl-1,3-benzothiazole"),
    ("Cc1cccc2scnc12",       "4-methyl-1,3-benzothiazole"),
    ("Cc1ccc2scnc2c1",       "5-methyl-1,3-benzothiazole"),
    # 1,3-benzoxazole: substitutable at C2, C4, C5, C6, C7
    ("c1ccc2ocnc2c1",        "1,3-benzoxazole"),
    ("Cc1nc2ccccc2o1",       "2-methyl-1,3-benzoxazole"),
    ("Cc1ccc2ocnc2c1",       "5-methyl-1,3-benzoxazole"),
    # 1,2-benzisothiazole: substitutable at C3, C4, C5, C6, C7
    ("c1ccc2sncc2c1",        "1,2-benzisothiazole"),
    ("Cc1ccc2sncc2c1",       "5-methyl-1,2-benzisothiazole"),
    # 1,2-benzisoxazole: substitutable at C3, C4, C5, C6, C7
    ("c1ccc2oncc2c1",        "1,2-benzisoxazole"),
    ("Cc1ccc2oncc2c1",       "5-methyl-1,2-benzisoxazole"),
    # 2,1-benzisothiazole: substitutable at C3, C4, C5, C6, C7
    ("c1ccc2nscc2c1",        "2,1-benzisothiazole"),
    ("Cc1ccc2nscc2c1",       "5-methyl-2,1-benzisothiazole"),
])
def test_phase554_benzo_azoles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
