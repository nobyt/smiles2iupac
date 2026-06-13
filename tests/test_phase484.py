"""Phase 484: imidazo/triazolo fused with [1,2,4]triazine and pyridazine
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cnn2cncc2n1",   "imidazo[1,5-b][1,2,4]triazine"),
    ("c1cnn2ccnc2n1",   "imidazo[3,2-b][1,2,4]triazine"),
    ("c1cnn2nncc2c1",   "[1,2,3]triazolo[1,5-b]pyridazine"),
])
def test_phase484_fused_triazine_pyridazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
