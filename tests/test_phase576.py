"""Phase 576: Substituted naphtho[2,1-b]furan and naphtho[2,1-b]thiophene naming.
Both share the same carbon skeleton (13 atoms); heteroatom at position 3.
Substitutable C positions: 1, 2, 4, 5, 6, 7, 8, 9.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # naphtho[2,1-b]furan
    ("c1ccc2c(c1)ccc1occc12",          "naphtho[2,1-b]furan"),
    ("Cc1coc2ccc3ccccc3c12",           "1-methylnaphtho[2,1-b]furan"),
    ("Cc1cc2c(ccc3ccccc32)o1",         "2-methylnaphtho[2,1-b]furan"),
    ("Cc1cc2ccccc2c2ccoc12",           "4-methylnaphtho[2,1-b]furan"),
    ("Cc1cc2occc2c2ccccc12",           "5-methylnaphtho[2,1-b]furan"),
    ("Cc1cccc2c1ccc1occc12",           "6-methylnaphtho[2,1-b]furan"),
    ("Cc1ccc2c(ccc3occc32)c1",         "7-methylnaphtho[2,1-b]furan"),
    ("Cc1ccc2ccc3occc3c2c1",           "8-methylnaphtho[2,1-b]furan"),
    ("Cc1cccc2ccc3occc3c12",           "9-methylnaphtho[2,1-b]furan"),
    # naphtho[2,1-b]thiophene
    ("c1ccc2c(c1)ccc1sccc12",          "naphtho[2,1-b]thiophene"),
    ("Cc1csc2ccc3ccccc3c12",           "1-methylnaphtho[2,1-b]thiophene"),
    ("Cc1cc2c(ccc3ccccc32)s1",         "2-methylnaphtho[2,1-b]thiophene"),
    ("Cc1cc2ccccc2c2ccsc12",           "4-methylnaphtho[2,1-b]thiophene"),
    ("Cc1cc2sccc2c2ccccc12",           "5-methylnaphtho[2,1-b]thiophene"),
    ("Cc1cccc2c1ccc1sccc12",           "6-methylnaphtho[2,1-b]thiophene"),
    ("Cc1ccc2c(ccc3sccc32)c1",         "7-methylnaphtho[2,1-b]thiophene"),
    ("Cc1ccc2ccc3sccc3c2c1",           "8-methylnaphtho[2,1-b]thiophene"),
    ("Cc1cccc2ccc3sccc3c12",           "9-methylnaphtho[2,1-b]thiophene"),
])
def test_phase576_naphtho_furan_thiophene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
