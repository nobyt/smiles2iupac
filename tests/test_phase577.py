"""Phase 577: Substituted naphtho[1,2-b]furan, naphtho[1,2-b]thiophene,
naphtho[2,3-b]furan, and naphtho[2,3-b]thiophene naming.
Heteroatom (O/S) is at position 1; substitutable C positions: 2-9.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # naphtho[1,2-b]furan
    ("c1ccc2c(c1)ccc1ccoc12",          "naphtho[1,2-b]furan"),
    ("Cc1cc2ccc3ccccc3c2o1",           "2-methylnaphtho[1,2-b]furan"),
    ("Cc1coc2c1ccc1ccccc12",           "3-methylnaphtho[1,2-b]furan"),
    ("Cc1cc2ccccc2c2occc12",           "4-methylnaphtho[1,2-b]furan"),
    ("Cc1cc2ccoc2c2ccccc12",           "5-methylnaphtho[1,2-b]furan"),
    ("Cc1cccc2c1ccc1ccoc12",           "6-methylnaphtho[1,2-b]furan"),
    ("Cc1ccc2c(ccc3ccoc32)c1",         "7-methylnaphtho[1,2-b]furan"),
    ("Cc1ccc2ccc3ccoc3c2c1",           "8-methylnaphtho[1,2-b]furan"),
    ("Cc1cccc2ccc3ccoc3c12",           "9-methylnaphtho[1,2-b]furan"),
    # naphtho[1,2-b]thiophene
    ("c1ccc2c(c1)ccc1ccsc12",          "naphtho[1,2-b]thiophene"),
    ("Cc1cc2ccc3ccccc3c2s1",           "2-methylnaphtho[1,2-b]thiophene"),
    ("Cc1csc2c1ccc1ccccc12",           "3-methylnaphtho[1,2-b]thiophene"),
    ("Cc1cc2ccccc2c2sccc12",           "4-methylnaphtho[1,2-b]thiophene"),
    ("Cc1cc2ccsc2c2ccccc12",           "5-methylnaphtho[1,2-b]thiophene"),
    ("Cc1cccc2c1ccc1ccsc12",           "6-methylnaphtho[1,2-b]thiophene"),
    ("Cc1ccc2c(ccc3ccsc32)c1",         "7-methylnaphtho[1,2-b]thiophene"),
    ("Cc1ccc2ccc3ccsc3c2c1",           "8-methylnaphtho[1,2-b]thiophene"),
    ("Cc1cccc2ccc3ccsc3c12",           "9-methylnaphtho[1,2-b]thiophene"),
    # naphtho[2,3-b]furan
    ("c1ccc2cc3occc3cc2c1",            "naphtho[2,3-b]furan"),
    ("Cc1cc2cc3ccccc3cc2o1",           "2-methylnaphtho[2,3-b]furan"),
    ("Cc1coc2cc3ccccc3cc12",           "3-methylnaphtho[2,3-b]furan"),
    ("Cc1c2ccccc2cc2occc12",           "4-methylnaphtho[2,3-b]furan"),
    ("Cc1cccc2cc3occc3cc12",           "5-methylnaphtho[2,3-b]furan"),
    ("Cc1ccc2cc3occc3cc2c1",           "6-methylnaphtho[2,3-b]furan"),
    ("Cc1ccc2cc3ccoc3cc2c1",           "7-methylnaphtho[2,3-b]furan"),
    ("Cc1cccc2cc3ccoc3cc12",           "8-methylnaphtho[2,3-b]furan"),
    ("Cc1c2ccccc2cc2ccoc12",           "9-methylnaphtho[2,3-b]furan"),
    # naphtho[2,3-b]thiophene
    ("c1ccc2cc3sccc3cc2c1",            "naphtho[2,3-b]thiophene"),
    ("Cc1cc2cc3ccccc3cc2s1",           "2-methylnaphtho[2,3-b]thiophene"),
    ("Cc1csc2cc3ccccc3cc12",           "3-methylnaphtho[2,3-b]thiophene"),
    ("Cc1c2ccccc2cc2sccc12",           "4-methylnaphtho[2,3-b]thiophene"),
    ("Cc1cccc2cc3sccc3cc12",           "5-methylnaphtho[2,3-b]thiophene"),
    ("Cc1ccc2cc3sccc3cc2c1",           "6-methylnaphtho[2,3-b]thiophene"),
    ("Cc1ccc2cc3ccsc3cc2c1",           "7-methylnaphtho[2,3-b]thiophene"),
    ("Cc1cccc2cc3ccsc3cc12",           "8-methylnaphtho[2,3-b]thiophene"),
    ("Cc1c2ccccc2cc2ccsc12",           "9-methylnaphtho[2,3-b]thiophene"),
])
def test_phase577_naphtho_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
