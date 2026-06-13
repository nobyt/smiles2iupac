"""Phase 500: remaining isothiazolo/isoxazolo orientations fused with [1,2,3]triazine at d-bond
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1nsc2cnnnc12",  "isothiazolo[4,5-d][1,2,3]triazine"),
    ("c1nnnc2oncc12",  "isoxazolo[5,4-d][1,2,3]triazine"),
])
def test_phase500_isothiazolo_isoxazolo_123triazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
