"""Phase 661: 3,4-dihydronaphthalen-1(2H)-one, indan-2-one methyl derivatives;
missing 3,4-dihydronaphthalen-2(1H)-one methyls (6,8)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 3,4-dihydronaphthalen-1(2H)-one (alpha-tetralone): C1(=O)-C2-C3-C4-C4a-...-C8a-C1
    # (parent in Phase 410; C1=O not methylable)
    ("O=C1CCCc2ccccc21",     "3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1C(C)CCc2ccccc21",  "2-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1CC(C)Cc2ccccc21",  "3-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1CCC(C)c2ccccc21",  "4-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1CCCc2c(C)cccc21",  "5-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1CCCc2cc(C)ccc21",  "6-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1CCCc2ccc(C)cc21",  "7-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1CCCc2cccc(C)c21",  "8-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    # 3,4-dihydronaphthalen-2(1H)-one (beta-tetralone): missing methyls 6,8
    # (Phase 628 covers positions 1,3,4,5,7; C2=O not methylable)
    ("O=C1CCc2cc(C)ccc2C1",  "6-methyl-3,4-dihydronaphthalen-2(1H)-one"),
    ("O=C1CCc2cccc(C)c2C1",  "8-methyl-3,4-dihydronaphthalen-2(1H)-one"),
    # indan-2-one: C1-C2(=O)-C3-C3a-C4-C5-C6-C7-C7a (C2-symmetric: 1≡3, 4≡7, 5≡6)
    # (parent in Phase 411; C2=O not methylable; lowest locants: 1 < 3, 4 < 7, 5 < 6)
    ("O=C1Cc2ccccc2C1",      "indan-2-one"),
    ("O=C1Cc2ccccc2C(C)1",   "1-methylindan-2-one"),
    ("Cc1cccc2c1CC(=O)C2",   "4-methylindan-2-one"),
    ("Cc1ccc2c(c1)CC(=O)C2", "5-methylindan-2-one"),
])
def test_phase661(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
