"""Phase 483: pyrazolo/triazolo fused with pyridazine and pyrazine
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cn2ncnc2cn1",   "[1,2,4]triazolo[1,5-a]pyrazine"),
    ("c1cnn2nccc2c1",   "pyrazolo[1,5-b]pyridazine"),
    ("c1cnn2ncnc2c1",   "[1,2,4]triazolo[1,5-b]pyridazine"),
])
def test_phase483_pyrazolo_triazolo_pyridazine_pyrazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
