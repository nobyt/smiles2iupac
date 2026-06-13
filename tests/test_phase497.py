"""Phase 497: 1H-pyrrolo and [1,2,3]triazolo fused with [1,2,3]triazine at d-bond
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2nnncc2[nH]1",  "1H-pyrrolo[3,2-d][1,2,3]triazine"),
    ("c1nnnc2c[nH]cc12",  "1H-pyrrolo[3,4-d][1,2,3]triazine"),
    ("c1cc2cnnnc2[nH]1",  "1H-pyrrolo[2,3-d][1,2,3]triazine"),
    ("c1nnnc2nn[nH]c12",  "1H-[1,2,3]triazolo[4,5-d][1,2,3]triazine"),
    ("c1nnnc2[nH]nnc12",  "1H-[1,2,3]triazolo[5,4-d][1,2,3]triazine"),
    ("c1nnnc2n[nH]nc12",  "2H-[1,2,3]triazolo[4,5-d][1,2,3]triazine"),
])
def test_phase497_pyrrolo_triazolo_123triazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
