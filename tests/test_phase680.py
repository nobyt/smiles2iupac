"""Phase 680: indan-2-one and benzazol-2(3H)-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # indan-2-one: C(1)-C(2,=O)-C(3)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 411; C2=O not methylable; C2v-symmetric: 1≡3, 4≡7, 5≡6)
    ("O=C1Cc2ccccc2C1",    "indan-2-one"),
    ("O=C1C(C)c2ccccc2C1", "1-methylindan-2-one"),
    ("O=C1Cc2c(C)cccc2C1", "4-methylindan-2-one"),
    ("O=C1Cc2cc(C)ccc2C1", "5-methylindan-2-one"),
    # 1H-benzimidazol-2(3H)-one: N(1,H)-C(2,=O)-N(3,H)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 412; C2=O not methylable; C2v-symmetric: 4≡7, 5≡6)
    ("O=c1[nH]c2ccccc2[nH]1",    "1H-benzimidazol-2(3H)-one"),
    ("O=c1[nH]c2c(C)cccc2[nH]1", "4-methyl-1H-benzimidazol-2(3H)-one"),
    ("O=c1[nH]c2cc(C)ccc2[nH]1", "5-methyl-1H-benzimidazol-2(3H)-one"),
    # 1,3-benzothiazol-2(3H)-one: S(1)-C(2,=O)-N(3,H)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 412; C2=O and S1 not methylable; N3-H and C4-C7 methylable)
    ("O=c1[nH]c2ccccc2s1",    "1,3-benzothiazol-2(3H)-one"),
    ("O=c1[nH]c2c(C)cccc2s1", "4-methyl-1,3-benzothiazol-2(3H)-one"),
    ("O=c1[nH]c2cc(C)ccc2s1", "5-methyl-1,3-benzothiazol-2(3H)-one"),
    ("O=c1[nH]c2ccc(C)cc2s1", "6-methyl-1,3-benzothiazol-2(3H)-one"),
    ("O=c1[nH]c2cccc(C)c2s1", "7-methyl-1,3-benzothiazol-2(3H)-one"),
    # 1,3-benzoxazol-2(3H)-one: O(1)-C(2,=O)-N(3,H)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 412; C2=O and O1 not methylable; N3-H and C4-C7 methylable)
    ("O=c1[nH]c2ccccc2o1",    "1,3-benzoxazol-2(3H)-one"),
    ("O=c1[nH]c2c(C)cccc2o1", "4-methyl-1,3-benzoxazol-2(3H)-one"),
    ("O=c1[nH]c2cc(C)ccc2o1", "5-methyl-1,3-benzoxazol-2(3H)-one"),
    ("O=c1[nH]c2ccc(C)cc2o1", "6-methyl-1,3-benzoxazol-2(3H)-one"),
    ("O=c1[nH]c2cccc(C)c2o1", "7-methyl-1,3-benzoxazol-2(3H)-one"),
])
def test_phase680(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
