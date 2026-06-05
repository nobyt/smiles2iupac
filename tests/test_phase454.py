"""Phase 454: benzo/pyrazino fused naphthyridines and pyrazinoquinoxalines
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cnc2cc3cccnc3cc2c1",    "benzo[f][1,5]naphthyridine"),
    ("c1ccc2c(c1)ncc1cnccc12",  "benzo[g][1,5]naphthyridine"),
    ("c1cnc2cc3nccnc3nc2c1",    "pyrazino[2,3-b][1,5]naphthyridine"),
    ("c1cnc2nc3nccnc3cc2c1",    "pyrazino[2,3-g][1,8]naphthyridine"),
    ("c1ccc2nc3nccnc3nc2c1",    "pyrazino[2,3-b]quinoxaline"),
    ("c1cnc2cc3nccnc3cc2n1",    "pyrazino[2,3-g]quinoxaline"),
    # regressions from Phase 453
    ("c1cc2cc3nccnc3cc2cn1",    "pyrido[3,4-g]quinoxaline"),
    ("c1cnc2cc3nccnc3cc2c1",    "pyrido[2,3-g]quinoxaline"),
    ("c1ccc2nc3nccnc3cc2c1",    "pyrido[2,3-b]quinoxaline"),
])
def test_phase454_fused_naphthyridines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
