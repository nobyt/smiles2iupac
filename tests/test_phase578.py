"""Phase 578: Substituted naphtho[2,3-d]oxazole, naphtho[2,3-d]thiazole,
naphtho[2,1-d]oxazole, and naphtho[2,1-d]thiazole naming.
Heteroatoms at positions 1 (O/S) and 3 (N); substitutable C positions: 2, 4-9.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # naphtho[2,3-d]oxazole
    ("c1ccc2cc3ocnc3cc2c1",            "naphtho[2,3-d]oxazole"),
    ("Cc1nc2cc3ccccc3cc2o1",           "2-methylnaphtho[2,3-d]oxazole"),
    ("Cc1c2ccccc2cc2ocnc12",           "4-methylnaphtho[2,3-d]oxazole"),
    ("Cc1cccc2cc3ocnc3cc12",           "5-methylnaphtho[2,3-d]oxazole"),
    ("Cc1ccc2cc3ocnc3cc2c1",           "6-methylnaphtho[2,3-d]oxazole"),
    ("Cc1ccc2cc3ncoc3cc2c1",           "7-methylnaphtho[2,3-d]oxazole"),
    ("Cc1cccc2cc3ncoc3cc12",           "8-methylnaphtho[2,3-d]oxazole"),
    ("Cc1c2ccccc2cc2ncoc12",           "9-methylnaphtho[2,3-d]oxazole"),
    # naphtho[2,3-d]thiazole
    ("c1ccc2cc3scnc3cc2c1",            "naphtho[2,3-d]thiazole"),
    ("Cc1nc2cc3ccccc3cc2s1",           "2-methylnaphtho[2,3-d]thiazole"),
    ("Cc1c2ccccc2cc2scnc12",           "4-methylnaphtho[2,3-d]thiazole"),
    ("Cc1cccc2cc3scnc3cc12",           "5-methylnaphtho[2,3-d]thiazole"),
    ("Cc1ccc2cc3scnc3cc2c1",           "6-methylnaphtho[2,3-d]thiazole"),
    ("Cc1ccc2cc3ncsc3cc2c1",           "7-methylnaphtho[2,3-d]thiazole"),
    ("Cc1cccc2cc3ncsc3cc12",           "8-methylnaphtho[2,3-d]thiazole"),
    ("Cc1c2ccccc2cc2ncsc12",           "9-methylnaphtho[2,3-d]thiazole"),
    # naphtho[2,1-d]oxazole
    ("c1ccc2c(c1)ccc1ncoc12",          "naphtho[2,1-d]oxazole"),
    ("Cc1nc2ccc3ccccc3c2o1",           "2-methylnaphtho[2,1-d]oxazole"),
    ("Cc1cc2ccccc2c2ocnc12",           "4-methylnaphtho[2,1-d]oxazole"),
    ("Cc1cc2ncoc2c2ccccc12",           "5-methylnaphtho[2,1-d]oxazole"),
    ("Cc1cccc2c1ccc1ncoc12",           "6-methylnaphtho[2,1-d]oxazole"),
    ("Cc1ccc2c(ccc3ncoc32)c1",         "7-methylnaphtho[2,1-d]oxazole"),
    ("Cc1ccc2ccc3ncoc3c2c1",           "8-methylnaphtho[2,1-d]oxazole"),
    ("Cc1cccc2ccc3ncoc3c12",           "9-methylnaphtho[2,1-d]oxazole"),
    # naphtho[2,1-d]thiazole
    ("c1ccc2c(c1)ccc1ncsc12",          "naphtho[2,1-d]thiazole"),
    ("Cc1nc2ccc3ccccc3c2s1",           "2-methylnaphtho[2,1-d]thiazole"),
    ("Cc1cc2ccccc2c2scnc12",           "4-methylnaphtho[2,1-d]thiazole"),
    ("Cc1cc2ncsc2c2ccccc12",           "5-methylnaphtho[2,1-d]thiazole"),
    ("Cc1cccc2c1ccc1ncsc12",           "6-methylnaphtho[2,1-d]thiazole"),
    ("Cc1ccc2c(ccc3ncsc32)c1",         "7-methylnaphtho[2,1-d]thiazole"),
    ("Cc1ccc2ccc3ncsc3c2c1",           "8-methylnaphtho[2,1-d]thiazole"),
    ("Cc1cccc2ccc3ncsc3c12",           "9-methylnaphtho[2,1-d]thiazole"),
])
def test_phase578_naphtho_oxazole_thiazole(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
