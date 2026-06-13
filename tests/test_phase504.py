"""Phase 504: [1,2,5]oxadiazolo fused bicyclics
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cnc2nonc2c1",   "[1,2,5]oxadiazolo[3,4-b]pyridine"),
    ("c1cnc2nonc2n1",   "[1,2,5]oxadiazolo[3,4-e]pyrazine"),
    ("c1nncc2nonc12",   "[1,2,5]oxadiazolo[3,4-d]pyridazine"),
    ("c1nnc2nonc2n1",   "[1,2,5]oxadiazolo[3,4-e][1,2,4]triazine"),
])
def test_phase504_125oxadiazolo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
