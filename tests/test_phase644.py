"""Phase 644: 3,4-dihydrophenanthridine, 1,2,3,4-tetrahydroacridine, 1,2,3,4-tetrahydrophenanthridine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 3,4-dihydrophenanthridine parent + methyl derivatives (OPSIN-verified; N5 aromatic, no 5-methyl)
    ("C1=Cc2c(ncc3ccccc23)CC1", "3,4-dihydrophenanthridine"),
    ("CC1=CCCc2ncc3ccccc3c21", "1-methyl-3,4-dihydrophenanthridine"),
    ("CC1=Cc2c(ncc3ccccc23)CC1", "2-methyl-3,4-dihydrophenanthridine"),
    ("CC1C=Cc2c(ncc3ccccc23)C1", "3-methyl-3,4-dihydrophenanthridine"),
    ("CC1CC=Cc2c1ncc1ccccc21", "4-methyl-3,4-dihydrophenanthridine"),
    ("Cc1nc2c(c3ccccc13)C=CCC2", "6-methyl-3,4-dihydrophenanthridine"),
    ("Cc1cccc2c3c(ncc12)CCC=C3", "7-methyl-3,4-dihydrophenanthridine"),
    ("Cc1ccc2c3c(ncc2c1)CCC=C3", "8-methyl-3,4-dihydrophenanthridine"),
    ("Cc1ccc2cnc3c(c2c1)C=CCC3", "9-methyl-3,4-dihydrophenanthridine"),
    ("Cc1cccc2cnc3c(c12)C=CCC3", "10-methyl-3,4-dihydrophenanthridine"),
    # 1,2,3,4-tetrahydroacridine parent + methyl derivatives (OPSIN-verified; N10 aromatic, no 10-methyl)
    ("c1ccc2nc3c(cc2c1)CCCC3", "1,2,3,4-tetrahydroacridine"),
    ("CC1CCCc2nc3ccccc3cc21", "1-methyl-1,2,3,4-tetrahydroacridine"),
    ("CC1CCc2nc3ccccc3cc2C1", "2-methyl-1,2,3,4-tetrahydroacridine"),
    ("CC1CCc2cc3ccccc3nc2C1", "3-methyl-1,2,3,4-tetrahydroacridine"),
    ("CC1CCCc2cc3ccccc3nc21", "4-methyl-1,2,3,4-tetrahydroacridine"),
    ("Cc1cccc2cc3c(nc12)CCCC3", "5-methyl-1,2,3,4-tetrahydroacridine"),
    ("Cc1ccc2cc3c(nc2c1)CCCC3", "6-methyl-1,2,3,4-tetrahydroacridine"),
    ("Cc1ccc2nc3c(cc2c1)CCCC3", "7-methyl-1,2,3,4-tetrahydroacridine"),
    ("Cc1cccc2nc3c(cc12)CCCC3", "8-methyl-1,2,3,4-tetrahydroacridine"),
    ("Cc1c2c(nc3ccccc13)CCCC2", "9-methyl-1,2,3,4-tetrahydroacridine"),
    # 1,2,3,4-tetrahydrophenanthridine parent + methyl derivatives (OPSIN-verified; N5 aromatic, no 5-methyl)
    ("c1ccc2c3c(ncc2c1)CCCC3", "1,2,3,4-tetrahydrophenanthridine"),
    ("CC1CCCc2ncc3ccccc3c21", "1-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("CC1CCc2ncc3ccccc3c2C1", "2-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("CC1CCc2c(ncc3ccccc23)C1", "3-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("CC1CCCc2c1ncc1ccccc21", "4-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("Cc1nc2c(c3ccccc13)CCCC2", "6-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("Cc1cccc2c3c(ncc12)CCCC3", "7-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("Cc1ccc2c3c(ncc2c1)CCCC3", "8-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("Cc1ccc2cnc3c(c2c1)CCCC3", "9-methyl-1,2,3,4-tetrahydrophenanthridine"),
    ("Cc1cccc2cnc3c(c12)CCCC3", "10-methyl-1,2,3,4-tetrahydrophenanthridine"),
])
def test_phase644(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
