"""Phase 493: 1H-imidazo and 1H-pyrazolo fused with [1,2,4]triazine at the e-bond
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1nnc2nc[nH]c2n1",  "1H-imidazo[5,4-e][1,2,4]triazine"),
    ("c1nnc2[nH]cnc2n1",  "1H-imidazo[4,5-e][1,2,4]triazine"),
    ("c1nnc2[nH]ncc2n1",  "1H-pyrazolo[4,5-e][1,2,4]triazine"),
    ("c1nnc2n[nH]cc2n1",  "1H-pyrazolo[4,3-e][1,2,4]triazine"),
    ("c1n[nH]c2cnnc-2n1",  "1H-pyrazolo[3,4-e][1,2,4]triazine"),
    ("c1n[nH]c2cnnc-2n1",  "1H-pyrazolo[3,4-e][1,2,4]triazine"),
])
def test_phase493_imidazo_pyrazolo_124triazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
