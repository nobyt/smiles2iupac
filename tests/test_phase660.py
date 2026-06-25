"""Phase 660: chroman-2-one, chroman-4-one, and isochroman-1-one methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # chroman-2-one (6H-1-benzopyran-2-one lactone): O1-C2(=O)-C3-C4-C4a-...-C8a-O1
    # (parent in Phase 410; C2=O and O1 not methylable)
    ("O=C1CCc2ccccc2O1",     "chroman-2-one"),
    ("CC1Cc2ccccc2OC1=O",    "3-methylchroman-2-one"),
    ("CC1CC(=O)Oc2ccccc21",  "4-methylchroman-2-one"),
    ("Cc1cccc2c1CCC(=O)O2",  "5-methylchroman-2-one"),
    ("Cc1ccc2c(c1)CCC(=O)O2","6-methylchroman-2-one"),
    ("Cc1ccc2c(c1)OC(=O)CC2","7-methylchroman-2-one"),
    ("Cc1cccc2c1OC(=O)CC2",  "8-methylchroman-2-one"),
    # chroman-4-one (4H-chromen-4-one precursor): O1-C2-C3-C4(=O)-C4a-...-C8a-O1
    # (parent in Phase 410; C4=O and O1 not methylable)
    ("O=C1CCOc2ccccc21",     "chroman-4-one"),
    ("O=C1CC(C)Oc2ccccc21",  "2-methylchroman-4-one"),
    ("O=C1C(C)COc2ccccc21",  "3-methylchroman-4-one"),
    ("O=C1CCOc2cccc(C)c21",  "5-methylchroman-4-one"),
    ("O=C1CCOc2ccc(C)cc21",  "6-methylchroman-4-one"),
    ("O=C1CCOc2cc(C)ccc21",  "7-methylchroman-4-one"),
    ("O=C1CCOc2c(C)cccc21",  "8-methylchroman-4-one"),
    # isochroman-1-one (3,4-dihydroisocoumarin): C1(=O)-O2-C3-C4-C4a-...-C8a-C1
    # (parent in Phase 410; C1=O and O2 not methylable)
    ("O=C1OCCc2ccccc21",     "isochroman-1-one"),
    ("O=C1OC(C)Cc2ccccc21",  "3-methylisochroman-1-one"),
    ("O=C1OCC(C)c2ccccc21",  "4-methylisochroman-1-one"),
    ("O=C1OCCc2c(C)cccc21",  "5-methylisochroman-1-one"),
    ("O=C1OCCc2cc(C)ccc21",  "6-methylisochroman-1-one"),
    ("O=C1OCCc2ccc(C)cc21",  "7-methylisochroman-1-one"),
    ("O=C1OCCc2cccc(C)c21",  "8-methylisochroman-1-one"),
])
def test_phase660(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
