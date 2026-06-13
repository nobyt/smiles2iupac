"""Phase 496: 1H-pyrrolo/[1,2,3]triazolo fused with [1,2,4]triazine at e-bond,
and 1H-pyrrolo[3,4-d]pyridazine (IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1nnc2cc[nH]c2n1",  "1H-pyrrolo[2,3-e][1,2,4]triazine"),
    ("c1nnc2c[nH]cc2n1",  "1H-pyrrolo[3,4-e][1,2,4]triazine"),
    ("c1nnc2[nH]ccc2n1",  "1H-pyrrolo[3,2-e][1,2,4]triazine"),
    ("c1nnc2nn[nH]c2n1",  "1H-[1,2,3]triazolo[5,4-e][1,2,4]triazine"),
    ("c1nnc2[nH]nnc2n1",  "1H-[1,2,3]triazolo[4,5-e][1,2,4]triazine"),
    ("c1nnc2n[nH]nc2n1",  "2H-[1,2,3]triazolo[4,5-e][1,2,4]triazine"),
    ("c1nncc2c[nH]cc12",  "1H-pyrrolo[3,4-d]pyridazine"),
])
def test_phase496_pyrrolo_triazolo_124triazine_and_pyridazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
