"""Phase 647: indoline and isoindoline (retained names) + methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # indoline parent + methyl derivatives (OPSIN-verified; locants 1-7)
    ("c1ccc2c(c1)CCN2", "2,3-dihydro-1H-indole"),
    ("CN1CCc2ccccc21", "1-methyl-2,3-dihydroindole"),
    ("CC1Cc2ccccc2N1", "2-methyl-2,3-dihydro-1H-indole"),
    ("CC1CNc2ccccc21", "3-methyl-2,3-dihydro-1H-indole"),
    ("Cc1cccc2c1CCN2", "4-methyl-2,3-dihydro-1H-indole"),
    ("Cc1ccc2c(c1)CCN2", "5-methyl-2,3-dihydro-1H-indole"),
    ("Cc1ccc2c(c1)NCC2", "6-methyl-2,3-dihydro-1H-indole"),
    ("Cc1cccc2c1NCC2", "7-methyl-2,3-dihydro-1H-indole"),
    # isoindoline parent + methyl derivatives (OPSIN-verified; C2-symmetric: 1=3, 4=7, 5=6)
    ("c1ccc2c(c1)CNC2", "2,3-dihydro-1H-isoindole"),
    ("CC1NCc2ccccc21", "1-methyl-2,3-dihydro-1H-isoindole"),
    ("CN1Cc2ccccc2C1", "2-methyl-1,3-dihydroisoindole"),
    ("Cc1cccc2c1CNC2", "4-methyl-2,3-dihydro-1H-isoindole"),
    ("Cc1ccc2c(c1)CNC2", "5-methyl-2,3-dihydro-1H-isoindole"),
])
def test_phase647(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
