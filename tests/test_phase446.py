"""Phase 446: benzo[f]quinoxaline, benzo[g]phthalazine, benzo[f]phthalazine
retained names (IUPAC 2013 P-31.1.3).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[f]quinoxaline (angular benzo+quinoxaline)
    ("c1ccc2c(c1)ccc1nccnc12",   "benzo[f]quinoxaline"),
    # benzo[g]phthalazine (linear benzo+phthalazine)
    ("c1ccc2cc3cnncc3cc2c1",     "benzo[g]phthalazine"),
    # benzo[f]phthalazine (angular benzo+phthalazine)
    ("c1ccc2c(c1)ccc1cnncc12",   "benzo[f]phthalazine"),
    # regressions
    ("c1ccc2cc3nccnc3cc2c1",     "benzo[g]quinoxaline"),
    ("c1ccc2cc3nnccc3cc2c1",     "benzo[g]cinnoline"),
    ("c1ccc2c(c1)ccc1ccnnc12",   "benzo[f]cinnoline"),
    ("c1ccc2nc3cnccc3cc2c1",     "benzo[c]cinnoline"),
    ("c1ccc2cnncc2c1",           "phthalazine"),
    ("c1ccc2nnccc2c1",           "cinnoline"),
])
def test_phase446_benzo_diaza(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
