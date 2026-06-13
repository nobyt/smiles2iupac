"""Phase 506: tetrazolo fused bicyclics
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cn2nnnc2cn1",   "tetrazolo[1,5-a]pyrazine"),
    ("c1cnn2nnnc2c1",   "tetrazolo[1,5-b]pyridazine"),
    ("c1cnn2nnnc2n1",   "tetrazolo[1,5-b][1,2,4]triazine"),
    ("c1cc2nnnn2nn1",   "tetrazolo[1,5-f][1,2,3]triazine"),
    ("c1nncn2nnnc12",   "tetrazolo[1,5-d][1,2,4]triazine"),
    ("c1ncc2nnnn2n1",   "tetrazolo[1,5-f][1,2,4]triazine"),
])
def test_phase506_tetrazolo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
