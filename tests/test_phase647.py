"""Phase 647: indoline and isoindoline (retained names) + methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # indoline parent + methyl derivatives (OPSIN-verified; locants 1-7)
    ("c1ccc2c(c1)CCN2", "indoline"),
    ("CN1CCc2ccccc21", "1-methylindoline"),
    ("CC1Cc2ccccc2N1", "2-methylindoline"),
    ("CC1CNc2ccccc21", "3-methylindoline"),
    ("Cc1cccc2c1CCN2", "4-methylindoline"),
    ("Cc1ccc2c(c1)CCN2", "5-methylindoline"),
    ("Cc1ccc2c(c1)NCC2", "6-methylindoline"),
    ("Cc1cccc2c1NCC2", "7-methylindoline"),
    # isoindoline parent + methyl derivatives (OPSIN-verified; C2-symmetric: 1=3, 4=7, 5=6)
    ("c1ccc2c(c1)CNC2", "isoindoline"),
    ("CC1NCc2ccccc21", "1-methylisoindoline"),
    ("CN1Cc2ccccc2C1", "2-methylisoindoline"),
    ("Cc1cccc2c1CNC2", "4-methylisoindoline"),
    ("Cc1ccc2c(c1)CNC2", "5-methylisoindoline"),
])
def test_phase647(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
