"""Phase 503: [1,2,3]thiadiazolo fused bicyclics
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cnc2nnsc2c1",   "[1,2,3]thiadiazolo[4,5-b]pyridine"),
    ("c1cnc2snnc2c1",   "[1,2,3]thiadiazolo[4,5-c]pyridine"),
    ("c1cnc2snnc2n1",   "[1,2,3]thiadiazolo[4,5-e]pyrazine"),
    ("c1nncc2snnc12",   "[1,2,3]thiadiazolo[5,4-d]pyridazine"),
    ("c1nnc2snnc2n1",   "[1,2,3]thiadiazolo[4,5-e][1,2,4]triazine"),
    ("c1nnnc2snnc12",   "[1,2,3]thiadiazolo[5,4-d][1,2,3]triazine"),
    ("c1nnnc2nnsc12",   "[1,2,3]thiadiazolo[4,5-d][1,2,3]triazine"),
])
def test_phase503_123thiadiazolo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
