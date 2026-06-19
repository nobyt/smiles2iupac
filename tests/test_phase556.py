"""Phase 556: Substituted 1,3-benzodioxole naming.
1,3-benzodioxole (methylenedioxy benzene) has locants 1(O), 2(CH2), 3(O),
3a and 7a (junctions), and 4–7 for the benzo ring carbons.
Positions 4-7 are substitutable; position 5 is most common (piperonal family).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # unsubstituted
    ("c1ccc2c(c1)OCO2",          "1,3-benzodioxole"),
    # C5-substituted (most common, piperonal family)
    ("Cc1ccc2c(c1)OCO2",         "5-methyl-1,3-benzodioxole"),
    ("Oc1ccc2c(c1)OCO2",         "1,3-benzodioxol-5-ol"),
    ("OC(=O)c1ccc2c(c1)OCO2",   "1,3-benzodioxole-5-carboxylic acid"),
    ("C=CCc1ccc2c(c1)OCO2",      "5-(prop-2-en-1-yl)-1,3-benzodioxole"),
    # C6-substituted
    ("c1c(C)cc2c(c1)OCO2",       "5-methyl-1,3-benzodioxole"),
    # C4-substituted (adjacent to junction)
    ("c1ccc2c(c1C)OCO2",         "4-methyl-1,3-benzodioxole"),
    # disubstituted
    ("Cc1cc2c(cc1C)OCO2",        "5,6-dimethyl-1,3-benzodioxole"),
])
def test_phase556_benzodioxole(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
