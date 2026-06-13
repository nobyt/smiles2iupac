"""Phase 489: isothiazolo/isoxazolo fused with [1,2,4]triazine at the e-bond
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1nnc2csnc2n1",  "isothiazolo[3,4-e][1,2,4]triazine"),
    ("c1nnc2nscc2n1",  "isothiazolo[4,3-e][1,2,4]triazine"),
    ("c1nnc2cnsc2n1",  "isothiazolo[5,4-e][1,2,4]triazine"),
    ("c1nnc2conc2n1",  "isoxazolo[3,4-e][1,2,4]triazine"),
    ("c1nnc2nocc2n1",  "isoxazolo[4,3-e][1,2,4]triazine"),
    ("c1nnc2cnoc2n1",  "isoxazolo[5,4-e][1,2,4]triazine"),
])
def test_phase489_isothiazolo_isoxazolo_triazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
