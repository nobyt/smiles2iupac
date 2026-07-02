"""Phase 746: 2-hydroxy-1,3-benzothiazole, 2-hydroxy-1,3-benzoxazole, and
2-hydroxy-1H-benzimidazole → their respective lactam tautomers (IUPAC 2013).

The enol (2-ol) forms prefer the N-H keto tautomers: 2(3H)-one.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Parent enol → keto tautomers
    ("Oc1nc2ccccc2s1",            "1,3-benzothiazol-2(3H)-one"),
    ("Oc1nc2ccccc2o1",            "1,3-benzoxazol-2(3H)-one"),
    ("Oc1nc2ccccc2[nH]1",         "1H-benzimidazol-2(3H)-one"),
    # Substituted: substituent prefix before lactam base
    ("Oc1nc2ccc(C)cc2s1",         "6-methyl-1,3-benzothiazol-2(3H)-one"),
    ("Oc1nc2ccc(C)cc2o1",         "6-methyl-1,3-benzoxazol-2(3H)-one"),
    # Regression: keto SMILES already correct
    ("O=c1[nH]c2ccccc2s1",        "1,3-benzothiazol-2(3H)-one"),
    ("O=c1[nH]c2ccccc2o1",        "1,3-benzoxazol-2(3H)-one"),
    ("O=c1[nH]c2ccccc2[nH]1",     "1H-benzimidazol-2(3H)-one"),
])
def test_phase746_benzoheterocycle_2ol_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
