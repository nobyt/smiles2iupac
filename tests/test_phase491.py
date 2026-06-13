"""Phase 491: thiazolo/oxazolo fused with [1,2,3]triazine at the d-bond
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1nc2cnnnc2s1",  "thiazolo[5,4-d][1,2,3]triazine"),
    ("c1nc2nnncc2s1",  "thiazolo[4,5-d][1,2,3]triazine"),
    ("c1nc2cnnnc2o1",  "oxazolo[5,4-d][1,2,3]triazine"),
    ("c1nc2nnncc2o1",  "oxazolo[4,5-d][1,2,3]triazine"),
])
def test_phase491_thiazolo_oxazolo_123triazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
