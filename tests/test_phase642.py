"""Phase 642: 1,2-dihydroacridine parent + methyl derivatives (locants 1-9)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Parent
    ("C1=Cc2nc3ccccc3cc2CC1", "1,2-dihydroacridine"),
    # Methyl derivatives (OPSIN-verified)
    ("CC1CC=Cc2nc3ccccc3cc21", "1-methyl-1,2-dihydroacridine"),
    ("CC1C=Cc2nc3ccccc3cc2C1", "2-methyl-1,2-dihydroacridine"),
    ("CC1=Cc2nc3ccccc3cc2CC1", "3-methyl-1,2-dihydroacridine"),
    ("CC1=CCCc2cc3ccccc3nc21", "4-methyl-1,2-dihydroacridine"),
    ("Cc1cccc2cc3c(nc12)C=CCC3", "5-methyl-1,2-dihydroacridine"),
    ("Cc1ccc2cc3c(nc2c1)C=CCC3", "6-methyl-1,2-dihydroacridine"),
    ("Cc1ccc2nc3c(cc2c1)CCC=C3", "7-methyl-1,2-dihydroacridine"),
    ("Cc1cccc2nc3c(cc12)CCC=C3", "8-methyl-1,2-dihydroacridine"),
    ("Cc1c2c(nc3ccccc13)C=CCC2", "9-methyl-1,2-dihydroacridine"),
])
def test_phase642(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
