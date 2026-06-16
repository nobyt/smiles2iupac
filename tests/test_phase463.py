"""Phase 463: 1,2,4-benzotriazine and pyrido[x,y-z]pyridazine isomers
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1ccc2nncnc2c1",   "1,2,4-benzotriazine"),
    ("c1cnc2ccnnc2c1",   "pyrido[2,3-e]pyridazine"),
    ("c1cnc2cnncc2c1",   "pyrido[2,3-d]pyridazine"),
    ("c1cc2ccnnc2cn1",   "pyrido[3,4-c]pyridazine"),
])
def test_phase463_benzotriazine_pyridopyridazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
