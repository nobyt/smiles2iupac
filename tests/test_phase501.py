"""Phase 501: missing isoxazolo[4,5-e][1,2,4]triazine
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1nnc2oncc2n1",  "isoxazolo[4,5-e][1,2,4]triazine"),
])
def test_phase501_isoxazolo_124triazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
