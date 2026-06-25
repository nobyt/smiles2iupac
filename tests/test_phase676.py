"""Phase 676: naphthalen-1(2H)-one, naphthalen-2(1H)-one, and 1,3-benzodioxol-2-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # naphthalen-1(2H)-one: C(1,=O)-C(2,sp3)-C(3)=C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent in Phase 419; C1=O not methylable; C2-C8 methylable)
    ("O=C1CC=Cc2ccccc21",    "naphthalen-1(2H)-one"),
    ("O=C1C(C)C=Cc2ccccc21", "2-methylnaphthalen-1(2H)-one"),
    ("O=C1CC(C)=Cc2ccccc21", "3-methylnaphthalen-1(2H)-one"),
    ("O=C1CC=C(C)c2ccccc21", "4-methylnaphthalen-1(2H)-one"),
    ("O=C1CC=Cc2c(C)cccc21", "5-methylnaphthalen-1(2H)-one"),
    ("O=C1CC=Cc2cc(C)ccc21", "6-methylnaphthalen-1(2H)-one"),
    ("O=C1CC=Cc2ccc(C)cc21", "7-methylnaphthalen-1(2H)-one"),
    ("O=C1CC=Cc2cccc(C)c21", "8-methylnaphthalen-1(2H)-one"),
    # naphthalen-2(1H)-one: C(2,=O)-C(3)=C(4)-C4a-C5-C6-C7-C8-C8a-C(1,sp3)
    # (parent in Phase 419; C2=O not methylable; C1 and C3-C8 methylable)
    ("O=C1C=Cc2ccccc2C1",    "naphthalen-2(1H)-one"),
    ("O=C1C=Cc2ccccc2C1C",   "1-methylnaphthalen-2(1H)-one"),
    ("O=C1C(C)=Cc2ccccc2C1", "3-methylnaphthalen-2(1H)-one"),
    ("O=C1C=C(C)c2ccccc2C1", "4-methylnaphthalen-2(1H)-one"),
    ("O=C1C=Cc2c(C)cccc2C1", "5-methylnaphthalen-2(1H)-one"),
    ("O=C1C=Cc2cc(C)ccc2C1", "6-methylnaphthalen-2(1H)-one"),
    ("O=C1C=Cc2ccc(C)cc2C1", "7-methylnaphthalen-2(1H)-one"),
    ("O=C1C=Cc2cccc(C)c2C1", "8-methylnaphthalen-2(1H)-one"),
    # 1,3-benzodioxol-2-one: O(1)-C(2,=O)-O(3)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 419; C2v-symmetric: 4≡7, 5≡6; all oxygens and C2 not methylable)
    ("O=c1oc2ccccc2o1",    "1,3-benzodioxol-2-one"),
    ("O=c1oc2c(C)cccc2o1", "4-methyl-1,3-benzodioxol-2-one"),
    ("O=c1oc2cc(C)ccc2o1", "5-methyl-1,3-benzodioxol-2-one"),
])
def test_phase676(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
