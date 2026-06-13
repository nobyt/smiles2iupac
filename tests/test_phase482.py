"""Phase 482: missing imidazo/pyrazolo/triazolo fused bicyclics
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2nccn2cn1",   "imidazo[1,2-c]pyrimidine"),
    ("c1cnc2cncn2c1",   "imidazo[1,5-a]pyrimidine"),
    ("c1cn2ncncc2n1",   "imidazo[1,2-f][1,2,4]triazine"),
    ("c1cn2nccc2cn1",   "pyrazolo[1,5-a]pyrazine"),
    ("c1cnn2cnnc2n1",   "[1,2,4]triazolo[4,3-b][1,2,4]triazine"),
])
def test_phase482_missing_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
