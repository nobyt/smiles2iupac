"""Phase 679: 3,4-dihydronaphthalen-1(2H)-one, isochroman-1-one, chroman-2-one, and chroman-4-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 3,4-dihydronaphthalen-1(2H)-one (beta-tetralone): C(1,=O)-C(2)-C(3)-C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent in Phase 410; C1=O not methylable; C2-C8 methylable)
    ("O=C1CCCc2ccccc21",    "3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1C(C)CCc2ccccc21", "2-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1CC(C)Cc2ccccc21", "3-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1CCC(C)c2ccccc21", "4-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1CCCc2c(C)cccc21", "5-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1CCCc2cc(C)ccc21", "6-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1CCCc2ccc(C)cc21", "7-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    ("O=C1CCCc2cccc(C)c21", "8-methyl-3,4-dihydronaphthalen-1(2H)-one"),
    # isochroman-1-one: C(1,=O)-O(2)-C(3)-C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent in Phase 410; C1=O and O2 not methylable; C3-C8 methylable)
    ("O=C1OCCc2ccccc21",    "isochroman-1-one"),
    ("O=C1OC(C)Cc2ccccc21", "3-methylisochroman-1-one"),
    ("O=C1OCC(C)c2ccccc21", "4-methylisochroman-1-one"),
    ("O=C1OCCc2c(C)cccc21", "5-methylisochroman-1-one"),
    ("O=C1OCCc2cc(C)ccc21", "6-methylisochroman-1-one"),
    ("O=C1OCCc2ccc(C)cc21", "7-methylisochroman-1-one"),
    ("O=C1OCCc2cccc(C)c21", "8-methylisochroman-1-one"),
    # chroman-2-one: O(1)-C(2,=O)-C(3)-C(4)-C4a-C5-C6-C7-C8-C8a
    # (parent in Phase 410; C2=O and O1 not methylable; C3-C8 methylable)
    ("O=C1CCc2ccccc2O1",    "chroman-2-one"),
    ("O=C1C(C)Cc2ccccc2O1", "3-methylchroman-2-one"),
    ("O=C1CC(C)c2ccccc2O1", "4-methylchroman-2-one"),
    ("O=C1CCc2c(C)cccc2O1", "5-methylchroman-2-one"),
    ("O=C1CCc2cc(C)ccc2O1", "6-methylchroman-2-one"),
    ("O=C1CCc2ccc(C)cc2O1", "7-methylchroman-2-one"),
    ("O=C1CCc2cccc(C)c2O1", "8-methylchroman-2-one"),
    # chroman-4-one: O(1)-C(2)-C(3)-C(4,=O)-C4a-C5-C6-C7-C8-C8a
    # (parent in Phase 410; C4=O and O1 not methylable; C2-C3 and C5-C8 methylable)
    ("O=C1CCOc2ccccc21",    "chroman-4-one"),
    ("O=C1CC(C)Oc2ccccc21", "2-methylchroman-4-one"),
    ("O=C1C(C)COc2ccccc21", "3-methylchroman-4-one"),
    ("O=C1CCOc2cccc(C)c21", "5-methylchroman-4-one"),
    ("O=C1CCOc2ccc(C)cc21", "6-methylchroman-4-one"),
    ("O=C1CCOc2cc(C)ccc21", "7-methylchroman-4-one"),
    ("O=C1CCOc2c(C)cccc21", "8-methylchroman-4-one"),
])
def test_phase679(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
