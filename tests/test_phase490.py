"""Phase 490: thieno/furo/isothiazolo/isoxazolo fused with [1,2,3]triazine at the d-bond
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2cnnnc2s1",  "thieno[2,3-d][1,2,3]triazine"),
    ("c1cc2nnncc2s1",  "thieno[3,2-d][1,2,3]triazine"),
    ("c1nnnc2cscc12",  "thieno[3,4-d][1,2,3]triazine"),
    ("c1cc2cnnnc2o1",  "furo[2,3-d][1,2,3]triazine"),
    ("c1cc2nnncc2o1",  "furo[3,2-d][1,2,3]triazine"),
    ("c1nnnc2cocc12",  "furo[3,4-d][1,2,3]triazine"),
    ("c1snc2cnnnc12",  "isothiazolo[4,3-d][1,2,3]triazine"),
    ("c1nnnc2nscc12",  "isothiazolo[3,4-d][1,2,3]triazine"),
    ("c1onc2cnnnc12",  "isoxazolo[4,3-d][1,2,3]triazine"),
    ("c1nnnc2nocc12",  "isoxazolo[3,4-d][1,2,3]triazine"),
])
def test_phase490_fused_123triazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
