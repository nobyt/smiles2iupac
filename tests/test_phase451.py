"""Phase 451: naphtho[1,2-d]oxazole/thiazole/pyrimidine and naphtho[1,2-b]furan/thiophene
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1ccc2c(c1)ccc1ocnc12",    "naphtho[1,2-d]oxazole"),
    ("c1ccc2c(c1)ccc1scnc12",    "naphtho[1,2-d]thiazole"),
    ("c1ccc2c(c1)ccc1ccoc12",    "naphtho[1,2-b]furan"),
    ("c1ccc2c(c1)ccc1ccsc12",    "naphtho[1,2-b]thiophene"),
    ("c1ccc2c(c1)ccc1cncnc12",   "naphtho[1,2-d]pyrimidine"),
    # regressions
    ("c1ccc2c(c1)ccc1[nH]cnc12", "1H-naphtho[2,1-d]imidazole"),
    ("c1ccc2c(c1)ccc1ncoc12",    "naphtho[2,1-d]oxazole"),
    ("c1ccc2c(c1)ccc1ncsc12",    "naphtho[2,1-d]thiazole"),
    ("c1ccc2c(c1)ccc1cn[nH]c12", "1H-naphtho[2,1-d]pyrazole"),
    ("c1ccc2cc3ncncc3cc2c1",     "naphtho[2,3-d]pyrimidine"),
])
def test_phase451_naphtho12d(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
