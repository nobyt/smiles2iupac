"""Phase 675: 1H-indole-2,3-dione, benzofuran-2(3H)-one, and benzo[b]thiophen-2(3H)-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-indole-2,3-dione (isatin): N(1,H)-C(2,=O)-C(3,=O)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 418; both C=O not methylable; N1-H methylable)
    ("O=C1Nc2ccccc2C1=O",    "1H-indole-2,3-dione"),
    ("O=C1N(C)c2ccccc2C1=O", "1-methylindole-2,3-dione"),
    ("O=C1Nc2c(C)cccc2C1=O", "4-methyl-1H-indole-2,3-dione"),
    ("O=C1Nc2cc(C)ccc2C1=O", "5-methyl-1H-indole-2,3-dione"),
    ("O=C1Nc2ccc(C)cc2C1=O", "6-methyl-1H-indole-2,3-dione"),
    ("O=C1Nc2cccc(C)c2C1=O", "7-methyl-1H-indole-2,3-dione"),
    # benzofuran-2(3H)-one: C(2,=O)-C(3,sp3)-O(1)-C7a-C4-C5-C6-C7-C3a
    # (parent in Phase 418; C2=O and O1 not methylable; C3-H2 and C4-C7 methylable)
    ("O=C1COc2ccccc21",    "1-benzofuran-2(3H)-one"),
    ("O=C1C(C)Oc2ccccc21", "3-methyl-1-benzofuran-2(3H)-one"),
    ("O=C1COc2cccc(C)c21", "4-methyl-1-benzofuran-2(3H)-one"),
    ("O=C1COc2ccc(C)cc21", "5-methyl-1-benzofuran-2(3H)-one"),
    ("O=C1COc2cc(C)ccc21", "6-methyl-1-benzofuran-2(3H)-one"),
    ("O=C1COc2c(C)cccc21", "7-methyl-1-benzofuran-2(3H)-one"),
    # benzo[b]thiophen-2(3H)-one: same layout with S at position 1
    # (parent in Phase 418; C2=O and S1 not methylable; C3-H2 and C4-C7 methylable)
    ("O=C1CSc2ccccc21",    "benzo[b]thiophen-2(3H)-one"),
    ("O=C1C(C)Sc2ccccc21", "3-methylbenzo[b]thiophen-2(3H)-one"),
    ("O=C1CSc2cccc(C)c21", "4-methylbenzo[b]thiophen-2(3H)-one"),
    ("O=C1CSc2ccc(C)cc21", "5-methylbenzo[b]thiophen-2(3H)-one"),
    ("O=C1CSc2cc(C)ccc21", "6-methylbenzo[b]thiophen-2(3H)-one"),
    ("O=C1CSc2c(C)cccc21", "7-methylbenzo[b]thiophen-2(3H)-one"),
])
def test_phase675(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
