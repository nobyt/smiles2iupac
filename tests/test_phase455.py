"""Phase 455: benzo/pyrido fused [1,6]-naphthyridines and pyrido-naphthyridines,
plus corrections from Phase 453 (IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1ccc2nc3ncccc3nc2c1",      "pyrido[2,3-b]quinoxaline"),
    ("c1cnc2c(c1)cnc1ccccc12",     "1,6-phenanthroline"),
    ("c1ccc2c(c1)ncc1ccncc12",    "2,6-phenanthroline"),
    ("c1cnc2nc3cccnc3cc2c1",      "pyrido[2,3-g][1,5]naphthyridine"),
    # regressions from Phase 453 fixes
    ("c1ccc2nc3nccnc3cc2c1",      "pyrazino[2,3-b]quinoline"),
    ("c1cnc2c(c1)cnc1nccnc12",      "pyrazino[2,3-h][1,6]naphthyridine"),
])
def test_phase455_fused_naphthyridines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
