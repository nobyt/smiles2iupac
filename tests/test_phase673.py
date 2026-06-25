"""Phase 673: indane-1,3-dione, isobenzofuran-1,3-dione, and fluoren-9-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # indane-1,3-dione: C(1,=O)-C(2)-C(3,=O)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 411; C2v-symmetric: 4≡7, 5≡6; both C=O not methylable)
    ("O=C1CC(=O)c2ccccc21",    "indane-1,3-dione"),
    ("O=C1C(C)C(=O)c2ccccc21", "2-methylindane-1,3-dione"),
    ("O=C1CC(=O)c2c(C)cccc21", "4-methylindane-1,3-dione"),
    ("O=C1CC(=O)c2cc(C)ccc21", "5-methylindane-1,3-dione"),
    # isobenzofuran-1,3-dione (phthalic anhydride): C(1,=O)-O(2)-C(3,=O)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 411; C2v-symmetric: 4≡7, 5≡6; both C=O and O not methylable)
    ("O=C1OC(=O)c2ccccc21",    "isobenzofuran-1,3-dione"),
    ("O=C1OC(=O)c2c(C)cccc21", "4-methylisobenzofuran-1,3-dione"),
    ("O=C1OC(=O)c2cc(C)ccc21", "5-methylisobenzofuran-1,3-dione"),
    # fluoren-9-one: C9(=O) / two benzo rings (C2v-symmetric: 1≡8, 2≡7, 3≡6, 4≡5)
    # (parent in Phase 408; C9=O not methylable; 4 unique positions)
    ("O=C1c2ccccc2-c2ccccc21",    "fluoren-9-one"),
    ("O=C1c2c(C)cccc2-c2ccccc21", "1-methylfluoren-9-one"),
    ("O=C1c2cc(C)ccc2-c2ccccc21", "2-methylfluoren-9-one"),
    ("O=C1c2ccc(C)cc2-c2ccccc21", "3-methylfluoren-9-one"),
    ("O=C1c2cccc(C)c2-c2ccccc21", "4-methylfluoren-9-one"),
])
def test_phase673(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
