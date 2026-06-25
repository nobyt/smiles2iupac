"""Phase 672: indan-1-one, indolin-2-one, and isoindolin-1-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # indan-1-one: C(1,=O)-C(2)-C(3)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 411; C1=O not methylable)
    ("O=C1CCc2ccccc21",    "indan-1-one"),
    ("O=C1C(C)Cc2ccccc21", "2-methylindan-1-one"),
    ("O=C1CC(C)c2ccccc21", "3-methylindan-1-one"),
    ("O=C1CCc2c(C)cccc21", "4-methylindan-1-one"),
    ("O=C1CCc2cc(C)ccc21", "5-methylindan-1-one"),
    ("O=C1CCc2ccc(C)cc21", "6-methylindan-1-one"),
    ("O=C1CCc2cccc(C)c21", "7-methylindan-1-one"),
    # indolin-2-one: N(1,H)-C(2,=O)-C(3)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 407; C2=O not methylable; N1-H methylable)
    ("O=C1Cc2ccccc2N1",    "indolin-2-one"),
    ("O=C1Cc2ccccc2N(C)1", "1-methylindolin-2-one"),
    ("O=C1C(C)c2ccccc2N1", "3-methylindolin-2-one"),
    ("O=C1Cc2c(C)cccc2N1", "4-methylindolin-2-one"),
    ("O=C1Cc2cc(C)ccc2N1", "5-methylindolin-2-one"),
    ("O=C1Cc2ccc(C)cc2N1", "6-methylindolin-2-one"),
    ("O=C1Cc2cccc(C)c2N1", "7-methylindolin-2-one"),
    # isoindolin-1-one: C(1,=O)-N(2,H)-C(3)-C3a-C4-C5-C6-C7-C7a
    # (parent in Phase 407; C1=O not methylable; N2-H methylable)
    ("O=C1NCc2ccccc21",    "isoindolin-1-one"),
    ("O=C1N(C)Cc2ccccc21", "2-methylisoindolin-1-one"),
    ("O=C1NC(C)c2ccccc21", "3-methylisoindolin-1-one"),
    ("O=C1NCc2c(C)cccc21", "4-methylisoindolin-1-one"),
    ("O=C1NCc2cc(C)ccc21", "5-methylisoindolin-1-one"),
    ("O=C1NCc2ccc(C)cc21", "6-methylisoindolin-1-one"),
    ("O=C1NCc2cccc(C)c21", "7-methylisoindolin-1-one"),
])
def test_phase672(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
