"""Phase 502: [1,2,3]oxadiazolo fused bicyclics
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cnc2nnoc2c1",   "[1,2,3]oxadiazolo[4,5-b]pyridine"),
    ("c1cc2onnc2cn1",   "[1,2,3]oxadiazolo[4,5-c]pyridine"),
    ("c1cnc2onnc2n1",   "[1,2,3]oxadiazolo[4,5-e]pyrazine"),
    ("c1nncc2onnc12",   "[1,2,3]oxadiazolo[5,4-d]pyridazine"),
    ("c1nnc2onnc2n1",   "[1,2,3]oxadiazolo[4,5-e][1,2,4]triazine"),
    ("c1nnnc2onnc12",   "[1,2,3]oxadiazolo[5,4-d][1,2,3]triazine"),
    ("c1nnnc2nnoc12",   "[1,2,3]oxadiazolo[4,5-d][1,2,3]triazine"),
])
def test_phase502_123oxadiazolo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
