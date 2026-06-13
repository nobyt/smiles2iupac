"""Phase 505: [1,2,5]thiadiazolo fused bicyclics
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cnc2nsnc2c1",   "[1,2,5]thiadiazolo[3,4-b]pyridine"),
    ("c1cnc2nsnc2n1",   "[1,2,5]thiadiazolo[3,4-e]pyrazine"),
    ("c1nncc2nsnc12",   "[1,2,5]thiadiazolo[3,4-d]pyridazine"),
    ("c1nnc2nsnc2n1",   "[1,2,5]thiadiazolo[3,4-e][1,2,4]triazine"),
])
def test_phase505_125thiadiazolo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
