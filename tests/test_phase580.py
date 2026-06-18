"""Phase 580: Substituted naphtho[2,3-d]isoxazole, naphtho[2,3-d]isothiazole,
naphtho[2,1-d]isoxazole, naphtho[2,1-d]isothiazole, naphtho[1,2-d]isoxazole,
and naphtho[1,2-d]isothiazole naming.
Heteroatoms: O/S at pos 1 (→ None), N at pos 2 (→ None); substitutable C: 3-9.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # naphtho[2,3-d]isoxazole
    ("c1ccc2cc3oncc3cc2c1",          "naphtho[2,3-d]isoxazole"),
    ("Cc1noc2cc3ccccc3cc12",         "3-methylnaphtho[2,3-d]isoxazole"),
    ("Cc1c2ccccc2cc2oncc12",         "4-methylnaphtho[2,3-d]isoxazole"),
    ("Cc1cccc2cc3oncc3cc12",         "5-methylnaphtho[2,3-d]isoxazole"),
    ("Cc1ccc2cc3oncc3cc2c1",         "6-methylnaphtho[2,3-d]isoxazole"),
    ("Cc1ccc2cc3cnoc3cc2c1",         "7-methylnaphtho[2,3-d]isoxazole"),
    ("Cc1cccc2cc3cnoc3cc12",         "8-methylnaphtho[2,3-d]isoxazole"),
    ("Cc1c2ccccc2cc2cnoc12",         "9-methylnaphtho[2,3-d]isoxazole"),
    # naphtho[2,3-d]isothiazole
    ("c1ccc2cc3sncc3cc2c1",          "naphtho[2,3-d]isothiazole"),
    ("Cc1nsc2cc3ccccc3cc12",         "3-methylnaphtho[2,3-d]isothiazole"),
    ("Cc1c2ccccc2cc2sncc12",         "4-methylnaphtho[2,3-d]isothiazole"),
    ("Cc1cccc2cc3sncc3cc12",         "5-methylnaphtho[2,3-d]isothiazole"),
    ("Cc1ccc2cc3sncc3cc2c1",         "6-methylnaphtho[2,3-d]isothiazole"),
    ("Cc1ccc2cc3cnsc3cc2c1",         "7-methylnaphtho[2,3-d]isothiazole"),
    ("Cc1cccc2cc3cnsc3cc12",         "8-methylnaphtho[2,3-d]isothiazole"),
    ("Cc1c2ccccc2cc2cnsc12",         "9-methylnaphtho[2,3-d]isothiazole"),
    # naphtho[2,1-d]isoxazole
    ("c1ccc2c(c1)ccc1cnoc12",        "naphtho[2,1-d]isoxazole"),
    ("Cc1noc2c1ccc1ccccc12",         "3-methylnaphtho[2,1-d]isoxazole"),
    ("Cc1cc2ccccc2c2oncc12",         "4-methylnaphtho[2,1-d]isoxazole"),
    ("Cc1cc2cnoc2c2ccccc12",         "5-methylnaphtho[2,1-d]isoxazole"),
    ("Cc1cccc2c1ccc1cnoc12",         "6-methylnaphtho[2,1-d]isoxazole"),
    ("Cc1ccc2c(ccc3cnoc32)c1",       "7-methylnaphtho[2,1-d]isoxazole"),
    ("Cc1ccc2ccc3cnoc3c2c1",         "8-methylnaphtho[2,1-d]isoxazole"),
    ("Cc1cccc2ccc3cnoc3c12",         "9-methylnaphtho[2,1-d]isoxazole"),
    # naphtho[2,1-d]isothiazole
    ("c1ccc2c(c1)ccc1cnsc12",        "naphtho[2,1-d]isothiazole"),
    ("Cc1nsc2c1ccc1ccccc12",         "3-methylnaphtho[2,1-d]isothiazole"),
    ("Cc1cc2ccccc2c2sncc12",         "4-methylnaphtho[2,1-d]isothiazole"),
    ("Cc1cc2cnsc2c2ccccc12",         "5-methylnaphtho[2,1-d]isothiazole"),
    ("Cc1cccc2c1ccc1cnsc12",         "6-methylnaphtho[2,1-d]isothiazole"),
    ("Cc1ccc2c(ccc3cnsc32)c1",       "7-methylnaphtho[2,1-d]isothiazole"),
    ("Cc1ccc2ccc3cnsc3c2c1",         "8-methylnaphtho[2,1-d]isothiazole"),
    ("Cc1cccc2ccc3cnsc3c12",         "9-methylnaphtho[2,1-d]isothiazole"),
    # naphtho[1,2-d]isoxazole
    ("c1ccc2c(c1)ccc1oncc12",        "naphtho[1,2-d]isoxazole"),
    ("Cc1noc2ccc3ccccc3c12",         "3-methylnaphtho[1,2-d]isoxazole"),
    ("Cc1cc2ccccc2c2cnoc12",         "4-methylnaphtho[1,2-d]isoxazole"),
    ("Cc1cc2oncc2c2ccccc12",         "5-methylnaphtho[1,2-d]isoxazole"),
    ("Cc1cccc2c1ccc1oncc12",         "6-methylnaphtho[1,2-d]isoxazole"),
    ("Cc1ccc2c(ccc3oncc32)c1",       "7-methylnaphtho[1,2-d]isoxazole"),
    ("Cc1ccc2ccc3oncc3c2c1",         "8-methylnaphtho[1,2-d]isoxazole"),
    ("Cc1cccc2ccc3oncc3c12",         "9-methylnaphtho[1,2-d]isoxazole"),
    # naphtho[1,2-d]isothiazole
    ("c1ccc2c(c1)ccc1sncc12",        "naphtho[1,2-d]isothiazole"),
    ("Cc1nsc2ccc3ccccc3c12",         "3-methylnaphtho[1,2-d]isothiazole"),
    ("Cc1cc2ccccc2c2cnsc12",         "4-methylnaphtho[1,2-d]isothiazole"),
    ("Cc1cc2sncc2c2ccccc12",         "5-methylnaphtho[1,2-d]isothiazole"),
    ("Cc1cccc2c1ccc1sncc12",         "6-methylnaphtho[1,2-d]isothiazole"),
    ("Cc1ccc2c(ccc3sncc32)c1",       "7-methylnaphtho[1,2-d]isothiazole"),
    ("Cc1ccc2ccc3sncc3c2c1",         "8-methylnaphtho[1,2-d]isothiazole"),
    ("Cc1cccc2ccc3sncc3c12",         "9-methylnaphtho[1,2-d]isothiazole"),
])
def test_phase580_naphtho_isoxazole_isothiazole(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
