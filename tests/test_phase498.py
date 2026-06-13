"""Phase 498: missing oxazolo/thiazolo/isothiazolo fused with [1,2,4]triazine at e-bond,
and isothiazolo/isoxazolo fused with [1,2,3]triazine at d-bond
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1nnc2sncc2n1",  "isothiazolo[4,5-e][1,2,4]triazine"),
    ("c1nnc2ocnc2n1",  "oxazolo[4,5-e][1,2,4]triazine"),
    ("c1nnc2ncoc2n1",  "oxazolo[5,4-e][1,2,4]triazine"),
    ("c1nnc2scnc2n1",  "thiazolo[4,5-e][1,2,4]triazine"),
    ("c1nnc2ncsc2n1",  "thiazolo[5,4-e][1,2,4]triazine"),
    ("c1nnnc2sncc12",  "isothiazolo[5,4-d][1,2,3]triazine"),
    ("c1noc2cnnnc12",  "isoxazolo[4,5-d][1,2,3]triazine"),
])
def test_phase498_oxazolo_thiazolo_triazines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
