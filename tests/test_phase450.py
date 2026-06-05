"""Phase 450: naphtho[2,1-d] series (imidazole, oxazole, thiazole, pyrazole)
and naphtho[2,3-d]pyrimidine (IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # naphtho[2,1-d] series
    ("c1ccc2c(c1)ccc1[nH]cnc12",  "1H-naphtho[2,1-d]imidazole"),
    ("c1ccc2c(c1)ccc1ocnc12",     "naphtho[2,1-d]oxazole"),
    ("c1ccc2c(c1)ccc1scnc12",     "naphtho[2,1-d]thiazole"),
    ("c1ccc2c(c1)ccc1[nH]ncc12",  "1H-naphtho[2,1-d]pyrazole"),
    # naphtho[2,3-d]pyrimidine
    ("c1ccc2cc3ncncc3cc2c1",      "naphtho[2,3-d]pyrimidine"),
    # regressions
    ("c1ccc2cc3[nH]cnc3cc2c1",    "1H-naphtho[2,3-d]imidazole"),
    ("c1ccc2cc3ocnc3cc2c1",       "naphtho[2,3-d]oxazole"),
    ("c1ccc2cc3scnc3cc2c1",       "naphtho[2,3-d]thiazole"),
    ("c1ccc2cc3[nH]ncc3cc2c1",    "1H-naphtho[2,3-d]pyrazole"),
    ("c1ccc2c(c1)nnc1ccccc12",    "benzo[c]cinnoline"),
])
def test_phase450_naphtho21d(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
