"""Phase 453: pyrido-quinoxaline and pyrido-quinoline tricyclics
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2cc3nccnc3cc2cn1",   "pyrido[3,4-g]quinoxaline"),
    ("c1cnc2cc3nccnc3cc2c1",   "pyrido[2,3-g]quinoxaline"),
    ("c1ccc2nc3nccnc3cc2c1",   "pyrazino[2,3-b]quinoline"),
    ("c1cc2nc3nccnc3cc2cn1",   "pyrazino[2,3-h][1,6]naphthyridine"),
    ("c1cnc2cc3ncccc3cc2c1",   "pyrido[2,3-g]quinoline"),
    # regressions
    ("c1ccc2cc3ncncc3cc2c1",   "naphtho[2,3-d]pyrimidine"),
    ("c1ccc2c(c1)ccc1ncncc12", "naphtho[1,2-d]pyrimidine"),
    ("c1ccc2c(c1)ccc1cncnc12", "naphtho[2,1-d]pyrimidine"),
    ("c1cnc2cc3ncccc3cc2c1",   "pyrido[2,3-g]quinoline"),
])
def test_phase453_pyrido_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
