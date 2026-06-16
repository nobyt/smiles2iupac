"""Phase 492: 1H-imidazo and 1H-pyrazolo fused with [1,2,3]triazine at the d-bond
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1nc2nnncc2[nH]1", "1H-imidazo[4,5-d][1,2,3]triazine"),
    ("c1nc2cnnnc2[nH]1", "1H-imidazo[5,4-d][1,2,3]triazine"),
    ("c1nnnc2[nH]ncc12", "1H-pyrazolo[5,4-d][1,2,3]triazine"),
    ("c1nnnc2n[nH]cc12", "1H-pyrazolo[3,4-d][1,2,3]triazine"),
    ("c1nn[nH]c2cnnc1-2", "1H-pyrazolo[4,3-d][1,2,3]triazine"),
    ("c1nn[nH]c2cnnc1-2", "1H-pyrazolo[4,3-d][1,2,3]triazine"),
])
def test_phase492_imidazo_pyrazolo_123triazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
