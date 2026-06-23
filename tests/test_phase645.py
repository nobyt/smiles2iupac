"""Phase 645: 1,2,3,4-tetrahydrobenzo[de]isoquinoline parent + methyl derivatives (locants 1-9)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Parent
    ("C1=Cc2cccc3c2C(C1)CNC3", "1,2,3,4-tetrahydrobenzo[de]isoquinoline"),
    # Methyl derivatives (OPSIN-verified; peri junction has no integer locant)
    ("CC1NCC2CC=Cc3cccc1c32", "1-methyl-1,2,3,4-tetrahydrobenzo[de]isoquinoline"),
    ("CN1Cc2cccc3c2C(CC=C3)C1", "2-methyl-1,2,3,4-tetrahydrobenzo[de]isoquinoline"),
    ("CC1NCc2cccc3c2C1CC=C3", "3-methyl-1,2,3,4-tetrahydrobenzo[de]isoquinoline"),
    ("CC1C=Cc2cccc3c2C1CNC3", "4-methyl-1,2,3,4-tetrahydrobenzo[de]isoquinoline"),
    ("CC1=Cc2cccc3c2C(CNC3)C1", "5-methyl-1,2,3,4-tetrahydrobenzo[de]isoquinoline"),
    ("CC1=CCC2CNCc3cccc1c32", "6-methyl-1,2,3,4-tetrahydrobenzo[de]isoquinoline"),
    ("Cc1ccc2c3c1C=CCC3CNC2", "7-methyl-1,2,3,4-tetrahydrobenzo[de]isoquinoline"),
    ("Cc1cc2c3c(c1)CNCC3CC=C2", "8-methyl-1,2,3,4-tetrahydrobenzo[de]isoquinoline"),
    ("Cc1ccc2c3c1CNCC3CC=C2", "9-methyl-1,2,3,4-tetrahydrobenzo[de]isoquinoline"),
])
def test_phase645(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
