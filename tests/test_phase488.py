"""Phase 488: thieno/furo fused with [1,2,4]triazine at the e-bond
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1nnc2ccsc2n1",  "thieno[2,3-e][1,2,4]triazine"),
    ("c1nnc2sccc2n1",  "thieno[3,2-e][1,2,4]triazine"),
    ("c1nnc2cscc2n1",  "thieno[3,4-e][1,2,4]triazine"),
    ("c1nnc2ccoc2n1",  "furo[2,3-e][1,2,4]triazine"),
    ("c1nnc2occc2n1",  "furo[3,2-e][1,2,4]triazine"),
    ("c1nnc2cocc2n1",  "furo[3,4-e][1,2,4]triazine"),
])
def test_phase488_thieno_furo_triazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
