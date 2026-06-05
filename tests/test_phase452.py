"""Phase 452: naphtho[2,3-b]furan/thiophene and naphtho[2,1-d]pyrimidine
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1ccc2cc3cocc3cc2c1",      "naphtho[2,3-b]furan"),
    ("c1ccc2cc3cscc3cc2c1",      "naphtho[2,3-b]thiophene"),
    ("c1ccc2c(c1)ccc1cncnc12",   "naphtho[2,1-d]pyrimidine"),
    # regressions
    ("c1ccc2c(c1)ccc1occc12",    "naphtho[1,2-b]furan"),
    ("c1ccc2c(c1)ccc1sccc12",    "naphtho[1,2-b]thiophene"),
    ("c1ccc2c(c1)ccc1ncncc12",   "naphtho[1,2-d]pyrimidine"),
    ("c1ccc2cc3ncncc3cc2c1",     "naphtho[2,3-d]pyrimidine"),
    ("c1ccc2c(c1)ccc1cnoc12",    "naphtho[1,2-d]oxazole"),
    ("c1ccc2c(c1)ccc1cnsc12",    "naphtho[1,2-d]thiazole"),
])
def test_phase452_naphtho23b(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
