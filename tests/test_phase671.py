"""Phase 671: isocoumarin and anthracen-9(10H)-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isocoumarin (1H-2-benzopyran-1-one): C(1,=O)-O(2)-C(3)-C(4)-C4a-...-C8a
    # (parent in Phase 550; C1=O and O2 not methylable)
    # note: ring traversal C1→O2→C3→C4 means first c-after-o = C3
    ("O=c1occc2ccccc12",    "isocoumarin"),
    ("O=c1oc(C)cc2ccccc12", "3-methylisocoumarin"),
    ("O=c1occ(C)c2ccccc12", "4-methylisocoumarin"),
    ("O=c1occc2c(C)cccc12", "5-methylisocoumarin"),
    ("O=c1occc2cc(C)ccc12", "6-methylisocoumarin"),
    ("O=c1occc2ccc(C)cc12", "7-methylisocoumarin"),
    ("O=c1occc2cccc(C)c12", "8-methylisocoumarin"),
    # anthracen-9(10H)-one (anthrone): C9(=O)-CH2(C10) / two benzo rings
    # (parent in Phase 415; C9=O and C10-H2 not methylable; C2v-symmetric: 1≡8, 2≡7, 3≡6, 4≡5)
    ("O=C1c2ccccc2Cc2ccccc21",    "anthracen-9(10H)-one"),
    ("O=C1c2c(C)cccc2Cc2ccccc21", "1-methylanthracen-9(10H)-one"),
    ("O=C1c2cc(C)ccc2Cc2ccccc21", "2-methylanthracen-9(10H)-one"),
    ("O=C1c2ccc(C)cc2Cc2ccccc21", "3-methylanthracen-9(10H)-one"),
    ("O=C1c2cccc(C)c2Cc2ccccc21", "4-methylanthracen-9(10H)-one"),
])
def test_phase671(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
