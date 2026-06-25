"""Phase 668: 1H-benzimidazol-2(3H)-one, 1,3-benzothiazol-2(3H)-one, and
1,3-benzoxazol-2(3H)-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-benzimidazol-2(3H)-one: N(1,H)-C(2,=O)-N(3,H)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 412; C2=O, N1-H, N3-H not methylable; C2v-symmetric: 4≡7, 5≡6)
    ("O=c1[nH]c2ccccc2[nH]1",    "1H-benzimidazol-2(3H)-one"),
    ("O=c1[nH]c2c(C)cccc2[nH]1", "4-methyl-1H-benzimidazol-2(3H)-one"),
    ("O=c1[nH]c2cc(C)ccc2[nH]1", "5-methyl-1H-benzimidazol-2(3H)-one"),
    # 1,3-benzothiazol-2(3H)-one: S(1)-C(2,=O)-N(3,H)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 412; C2=O and N3-H not methylable)
    ("O=c1[nH]c2ccccc2s1",        "1,3-benzothiazol-2(3H)-one"),
    ("O=c1[nH]c2c(C)cccc2s1",     "4-methyl-1,3-benzothiazol-2(3H)-one"),
    ("O=c1[nH]c2cc(C)ccc2s1",     "5-methyl-1,3-benzothiazol-2(3H)-one"),
    ("O=c1[nH]c2ccc(C)cc2s1",     "6-methyl-1,3-benzothiazol-2(3H)-one"),
    ("O=c1[nH]c2cccc(C)c2s1",     "7-methyl-1,3-benzothiazol-2(3H)-one"),
    # 1,3-benzoxazol-2(3H)-one: O(1)-C(2,=O)-N(3,H)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 412; C2=O and N3-H not methylable)
    ("O=c1[nH]c2ccccc2o1",        "1,3-benzoxazol-2(3H)-one"),
    ("O=c1[nH]c2c(C)cccc2o1",     "4-methyl-1,3-benzoxazol-2(3H)-one"),
    ("O=c1[nH]c2cc(C)ccc2o1",     "5-methyl-1,3-benzoxazol-2(3H)-one"),
    ("O=c1[nH]c2ccc(C)cc2o1",     "6-methyl-1,3-benzoxazol-2(3H)-one"),
    ("O=c1[nH]c2cccc(C)c2o1",     "7-methyl-1,3-benzoxazol-2(3H)-one"),
])
def test_phase668(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
