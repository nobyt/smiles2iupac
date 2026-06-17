"""Phase 565: Substituted benz[a]acridine and benz[c]acridine naming.
Both are 18-atom tetracyclic 6+6+6+6 linearly fused ring systems with N in one ring.
benz[a]acridine: N at position 10; substitutable C positions 1-9,11,12.
benz[c]acridine: N at position 7; substitutable C positions 1-6,8-12.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benz[a]acridine (N at pos 10; C positions 1-9,11,12)
    ("c1ccc2nc3ccc4ccccc4c3cc2c1",          "benz[a]acridine"),
    ("Cc1cccc2c1ccc1nc3ccccc3cc12",          "1-methylbenz[a]acridine"),
    ("Cc1ccc2c(ccc3nc4ccccc4cc32)c1",        "2-methylbenz[a]acridine"),
    ("Cc1ccc2ccc3nc4ccccc4cc3c2c1",          "3-methylbenz[a]acridine"),
    ("Cc1cccc2ccc3nc4ccccc4cc3c12",          "4-methylbenz[a]acridine"),
    ("Cc1c2ccccc2nc2ccc3ccccc3c12",          "5-methylbenz[a]acridine"),
    ("Cc1cccc2nc3ccc4ccccc4c3cc12",          "6-methylbenz[a]acridine"),
    ("Cc1ccc2nc3ccc4ccccc4c3cc2c1",          "7-methylbenz[a]acridine"),
    ("Cc1ccc2cc3c(ccc4ccccc43)nc2c1",        "8-methylbenz[a]acridine"),
    ("Cc1cccc2cc3c(ccc4ccccc43)nc12",        "9-methylbenz[a]acridine"),
    ("Cc1cc2ccccc2c2cc3ccccc3nc12",          "11-methylbenz[a]acridine"),
    ("Cc1cc2nc3ccccc3cc2c2ccccc12",          "12-methylbenz[a]acridine"),
    # benz[c]acridine (N at pos 7; C positions 1-6,8-12)
    ("c1ccc2cc3nc4ccccc4cc3cc2c1",           "benz[c]acridine"),
    ("Cc1c2ccccc2cc2nc3ccccc3cc12",          "1-methylbenz[c]acridine"),
    ("Cc1cccc2cc3nc4ccccc4cc3cc12",          "2-methylbenz[c]acridine"),
    ("Cc1ccc2cc3nc4ccccc4cc3cc2c1",          "3-methylbenz[c]acridine"),
    ("Cc1ccc2cc3cc4ccccc4nc3cc2c1",          "4-methylbenz[c]acridine"),
    ("Cc1cccc2cc3cc4ccccc4nc3cc12",          "5-methylbenz[c]acridine"),
    ("Cc1c2ccccc2cc2cc3ccccc3nc12",          "6-methylbenz[c]acridine"),
    ("Cc1cccc2cc3cc4ccccc4cc3nc12",          "8-methylbenz[c]acridine"),
    ("Cc1ccc2cc3cc4ccccc4cc3nc2c1",          "9-methylbenz[c]acridine"),
    ("Cc1ccc2nc3cc4ccccc4cc3cc2c1",          "10-methylbenz[c]acridine"),
    ("Cc1cccc2nc3cc4ccccc4cc3cc12",          "11-methylbenz[c]acridine"),
    ("Cc1c2ccccc2nc2cc3ccccc3cc12",          "12-methylbenz[c]acridine"),
])
def test_phase565_benz_acridines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
