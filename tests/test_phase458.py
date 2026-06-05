"""Phase 458: pyrimido[4,5-b]quinoxaline, pyrido[2,3-b][1,6/7/8]naphthyridines
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrimidine at bond b of quinoxaline; N3 adjacent to C2 (lower locant)
    ("c1ccc2nc3ncncc3nc2c1",      "pyrimido[4,5-b]quinoxaline"),
    # pyridine at bond b of [1,8]naphthyridine; N adjacent to C2 (lower locant)
    ("c1cnc2nc3ncccc3cc2c1",      "pyrido[2,3-b][1,8]naphthyridine"),
    # pyridine at bond b of [1,6]naphthyridine; N adjacent to C2 (lower locant)
    ("c1cnc2nc3ccncc3cc2c1",      "pyrido[2,3-b][1,6]naphthyridine"),
    # pyridine at bond b of [1,7]naphthyridine; N adjacent to C2 (lower locant)
    ("c1cnc2nc3cnccc3cc2c1",      "pyrido[2,3-b][1,7]naphthyridine"),
])
def test_phase458_pyrido_naphthyridines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
