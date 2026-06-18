"""Phase 583: Substituted naphtho[2,3-d]pyrimidine, naphtho[2,1-d]pyrimidine,
and naphtho[1,2-d]pyrimidine naming.
N atoms at positions 1&3 (2,1 isomer) or 1&3 at different ring positions.
naphtho[2,3-d]: sub C 2,4,5,6,7,8,9,10; naphtho[2,1-d]: sub C 1,3,5,6,7,8,9,10;
naphtho[1,2-d]: sub C 2,4,5,6,7,8,9,10.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # naphtho[2,3-d]pyrimidine
    ("c1ccc2cc3ncncc3cc2c1",         "naphtho[2,3-d]pyrimidine"),
    ("Cc1ncc2cc3ccccc3cc2n1",        "2-methylnaphtho[2,3-d]pyrimidine"),
    ("Cc1ncnc2cc3ccccc3cc12",        "4-methylnaphtho[2,3-d]pyrimidine"),
    ("Cc1c2ccccc2cc2ncncc12",        "5-methylnaphtho[2,3-d]pyrimidine"),
    ("Cc1cccc2cc3ncncc3cc12",        "6-methylnaphtho[2,3-d]pyrimidine"),
    ("Cc1ccc2cc3ncncc3cc2c1",        "7-methylnaphtho[2,3-d]pyrimidine"),
    ("Cc1ccc2cc3cncnc3cc2c1",        "8-methylnaphtho[2,3-d]pyrimidine"),
    ("Cc1cccc2cc3cncnc3cc12",        "9-methylnaphtho[2,3-d]pyrimidine"),
    ("Cc1c2ccccc2cc2cncnc12",        "10-methylnaphtho[2,3-d]pyrimidine"),
    # naphtho[2,1-d]pyrimidine
    ("c1ccc2c(c1)ccc1ncncc12",       "naphtho[2,1-d]pyrimidine"),
    ("Cc1ncnc2ccc3ccccc3c12",        "1-methylnaphtho[2,1-d]pyrimidine"),
    ("Cc1ncc2c(ccc3ccccc32)n1",      "3-methylnaphtho[2,1-d]pyrimidine"),
    ("Cc1cc2ccccc2c2cncnc12",        "5-methylnaphtho[2,1-d]pyrimidine"),
    ("Cc1cc2ncncc2c2ccccc12",        "6-methylnaphtho[2,1-d]pyrimidine"),
    ("Cc1cccc2c1ccc1ncncc12",        "7-methylnaphtho[2,1-d]pyrimidine"),
    ("Cc1ccc2c(ccc3ncncc32)c1",      "8-methylnaphtho[2,1-d]pyrimidine"),
    ("Cc1ccc2ccc3ncncc3c2c1",        "9-methylnaphtho[2,1-d]pyrimidine"),
    ("Cc1cccc2ccc3ncncc3c12",        "10-methylnaphtho[2,1-d]pyrimidine"),
    # naphtho[1,2-d]pyrimidine
    ("c1ccc2c(c1)ccc1cncnc12",       "naphtho[1,2-d]pyrimidine"),
    ("Cc1ncc2ccc3ccccc3c2n1",        "2-methylnaphtho[1,2-d]pyrimidine"),
    ("Cc1ncnc2c1ccc1ccccc12",        "4-methylnaphtho[1,2-d]pyrimidine"),
    ("Cc1cc2ccccc2c2ncncc12",        "5-methylnaphtho[1,2-d]pyrimidine"),
    ("Cc1cc2cncnc2c2ccccc12",        "6-methylnaphtho[1,2-d]pyrimidine"),
    ("Cc1cccc2c1ccc1cncnc12",        "7-methylnaphtho[1,2-d]pyrimidine"),
    ("Cc1ccc2c(ccc3cncnc32)c1",      "8-methylnaphtho[1,2-d]pyrimidine"),
    ("Cc1ccc2ccc3cncnc3c2c1",        "9-methylnaphtho[1,2-d]pyrimidine"),
    ("Cc1cccc2ccc3cncnc3c12",        "10-methylnaphtho[1,2-d]pyrimidine"),
])
def test_phase583_naphtho_pyrimidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
