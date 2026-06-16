"""Phase 553: Substituted 1H-indazole and 1H-benzotriazole naming.
1H-indazole (benzo[c]pyrazole) has locants 1(NH),2(N),3,3a,4,5,6,7,7a.
1H-benzotriazole has locants 1(NH),2(N),3(N),3a,4,5,6,7,7a.
Both share the same 9-atom locant map as indole/benzimidazole/benzofuran.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-indazole unsubstituted
    ("c1ccc2[nH]ncc2c1",         "1H-indazole"),
    # substituted 1H-indazole
    ("Cc1n[nH]c2ccccc12",        "3-methyl-1H-indazole"),
    ("Cc1cccc2[nH]ncc12",        "4-methyl-1H-indazole"),
    ("Cc1ccc2[nH]ncc2c1",        "5-methyl-1H-indazole"),
    ("c1c(C)cc2[nH]ncc2c1",      "6-methyl-1H-indazole"),
    ("Cc1cccc2c1[nH]nc2",        "7-methyl-1H-indazole"),
    # 1H-benzotriazole unsubstituted
    ("c1ccc2[nH]nnc2c1",         "1H-benzotriazole"),
    # substituted 1H-benzotriazole
    ("Cc1ccc2[nH]nnc2c1",        "5-methyl-1H-benzotriazole"),
    ("Cc1cccc2[nH]nnc12",        "4-methyl-1H-benzotriazole"),
])
def test_phase553_indazole_benzotriazole(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
