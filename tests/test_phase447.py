"""Phase 447: benzo[h]cinnoline, furo[3,2-b]quinoline, furo[3,4-b]quinoline
retained names (IUPAC 2013 P-31.1.3).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[h]cinnoline
    ("c1ccc2c(c1)ccc1nnccc12",   "benzo[h]cinnoline"),
    # furo[3,2-b]quinoline (O adjacent to quinoline-C3)
    ("c1ccc2nc3ccoc3cc2c1",      "furo[3,2-b]quinoline"),
    # furo[3,4-b]quinoline
    ("c1ccc2nc3cocc3cc2c1",      "furo[3,4-b]quinoline"),
    # regressions
    ("c1ccc2c(c1)ccc1ccnnc12",   "benzo[f]cinnoline"),
    ("c1ccc2cc3nnccc3cc2c1",     "benzo[g]cinnoline"),
    ("c1ccc2nc3cnccc3cc2c1",     "benzo[b][1,7]naphthyridine"),
    ("c1ccc2c(c1)ccc1nccnc12",   "benzo[f]quinoxaline"),
    ("c1ccc2cc3nccnc3cc2c1",     "benzo[g]quinoxaline"),
])
def test_phase447_benzo_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
