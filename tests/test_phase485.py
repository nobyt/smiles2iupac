"""Phase 485: pyrrolo/imidazo/triazolo/pyrazolo fused with pyrazine/pyridazine
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2cnccn2c1",  "pyrrolo[1,2-a]pyrazine"),
    ("c1cnn2cccc2c1",  "pyrrolo[1,2-b]pyridazine"),
    ("c1cn2cncc2cn1",  "imidazo[1,5-a]pyrazine"),
    ("c1cnn2cncc2c1",  "imidazo[1,5-b]pyridazine"),
    ("c1cn2cnnc2cn1",  "[1,2,4]triazolo[4,3-a]pyrazine"),
    ("c1cnn2nccc2n1",  "pyrazolo[1,5-b][1,2,4]triazine"),
])
def test_phase485_fused_pyrazine_pyridazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
