"""Phase 464: pyrido[x,y-e]pyridazine, pyrido[x,y-e]pyrimidine, pyrido[x,y-e]pyrazine
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2nnccc2cn1",   "pyrido[3,4-e]pyridazine"),
    ("c1cnc2nnccc2c1",   "pyrido[2,3-e]pyridazine"),
    ("c1cnc2ncncc2c1",   "pyrido[2,3-e]pyrimidine"),
    ("c1cc2ncncc2cn1",   "pyrido[3,4-e]pyrimidine"),
    ("c1cc2nccnc2cn1",   "pyrido[3,4-e]pyrazine"),
])
def test_phase464_pyrido_diazine_isomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
