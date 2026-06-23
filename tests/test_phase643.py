"""Phase 643: 3,4-dihydroacridine (locants 1-9) and 1,2-dihydrophenanthridine (locants 1-4, 6-10)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 3,4-dihydroacridine parent + methyl derivatives (OPSIN-verified)
    ("C1=Cc2cc3ccccc3nc2CC1", "3,4-dihydroacridine"),
    ("CC1=CCCc2nc3ccccc3cc21", "1-methyl-3,4-dihydroacridine"),
    ("CC1=Cc2cc3ccccc3nc2CC1", "2-methyl-3,4-dihydroacridine"),
    ("CC1C=Cc2cc3ccccc3nc2C1", "3-methyl-3,4-dihydroacridine"),
    ("CC1CC=Cc2cc3ccccc3nc21", "4-methyl-3,4-dihydroacridine"),
    ("Cc1cccc2cc3c(nc12)CCC=C3", "5-methyl-3,4-dihydroacridine"),
    ("Cc1ccc2cc3c(nc2c1)CCC=C3", "6-methyl-3,4-dihydroacridine"),
    ("Cc1ccc2nc3c(cc2c1)C=CCC3", "7-methyl-3,4-dihydroacridine"),
    ("Cc1cccc2nc3c(cc12)C=CCC3", "8-methyl-3,4-dihydroacridine"),
    ("Cc1c2c(nc3ccccc13)CCC=C2", "9-methyl-3,4-dihydroacridine"),
    # 1,2-dihydrophenanthridine parent + methyl derivatives (OPSIN-verified; N5 aromatic, no 5-methyl)
    ("C1=Cc2ncc3ccccc3c2CC1", "1,2-dihydrophenanthridine"),
    ("CC1CC=Cc2ncc3ccccc3c21", "1-methyl-1,2-dihydrophenanthridine"),
    ("CC1C=Cc2ncc3ccccc3c2C1", "2-methyl-1,2-dihydrophenanthridine"),
    ("CC1=Cc2ncc3ccccc3c2CC1", "3-methyl-1,2-dihydrophenanthridine"),
    ("CC1=CCCc2c1ncc1ccccc21", "4-methyl-1,2-dihydrophenanthridine"),
    ("Cc1nc2c(c3ccccc13)CCC=C2", "6-methyl-1,2-dihydrophenanthridine"),
    ("Cc1cccc2c3c(ncc12)C=CCC3", "7-methyl-1,2-dihydrophenanthridine"),
    ("Cc1ccc2c3c(ncc2c1)C=CCC3", "8-methyl-1,2-dihydrophenanthridine"),
    ("Cc1ccc2cnc3c(c2c1)CCC=C3", "9-methyl-1,2-dihydrophenanthridine"),
    ("Cc1cccc2cnc3c(c12)CCC=C3", "10-methyl-1,2-dihydrophenanthridine"),
])
def test_phase643(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
