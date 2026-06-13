"""Phase 487: thieno/furo pyridazine [3,4-c] and [x,y-d] isomers
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2cscc2nn1",  "thieno[3,4-c]pyridazine"),
    ("c1cc2cocc2nn1",  "furo[3,4-c]pyridazine"),
    ("c1cc2cnncc2s1",  "thieno[3,2-d]pyridazine"),
    ("c1nncc2cscc12",  "thieno[3,4-d]pyridazine"),
    ("c1cc2cnncc2o1",  "furo[3,2-d]pyridazine"),
    ("c1nncc2cocc12",  "furo[3,4-d]pyridazine"),
])
def test_phase487_thieno_furo_pyridazine_isomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
