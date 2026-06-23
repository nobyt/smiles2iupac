"""Phase 646: 5,6,7,8-tetrahydrophenanthridine parent + methyl derivatives (locants 1-10)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Parent
    ("C1=CC2=C(CC1)CNc1ccccc12", "5,6,7,8-tetrahydrophenanthridine"),
    # Methyl derivatives (OPSIN-verified)
    ("Cc1cccc2c1C1=C(CCC=C1)CN2", "1-methyl-5,6,7,8-tetrahydrophenanthridine"),
    ("Cc1ccc2c(c1)C1=C(CCC=C1)CN2", "2-methyl-5,6,7,8-tetrahydrophenanthridine"),
    ("Cc1ccc2c(c1)NCC1=C2C=CCC1", "3-methyl-5,6,7,8-tetrahydrophenanthridine"),
    ("Cc1cccc2c1NCC1=C2C=CCC1", "4-methyl-5,6,7,8-tetrahydrophenanthridine"),
    ("CN1CC2=C(C=CCC2)c2ccccc21", "5-methyl-5,6,7,8-tetrahydrophenanthridine"),
    ("CC1Nc2ccccc2C2=C1CCC=C2", "6-methyl-5,6,7,8-tetrahydrophenanthridine"),
    ("CC1CC=CC2=C1CNc1ccccc12", "7-methyl-5,6,7,8-tetrahydrophenanthridine"),
    ("CC1C=CC2=C(CNc3ccccc32)C1", "8-methyl-5,6,7,8-tetrahydrophenanthridine"),
    ("CC1=CC2=C(CC1)CNc1ccccc12", "9-methyl-5,6,7,8-tetrahydrophenanthridine"),
    ("CC1=CCCC2=C1c1ccccc1NC2", "10-methyl-5,6,7,8-tetrahydrophenanthridine"),
])
def test_phase646(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
