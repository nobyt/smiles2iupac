"""Phase 494: 1H-pyrrolo/pyrazolo/imidazo fused with pyridazine at c- and d-bonds
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2cn[nH]c2nn1",  "1H-pyrazolo[5,4-c]pyridazine"),
    ("c1cc2n[nH]cc2nn1",  "1H-pyrazolo[4,3-c]pyridazine"),
    ("c1nc2cnncc2[nH]1",  "1H-imidazo[4,5-d]pyridazine"),
    ("c1nncc2[nH]ncc12",  "1H-pyrazolo[4,5-d]pyridazine"),
    ("c1nncc2n[nH]cc12",  "1H-pyrazolo[3,4-d]pyridazine"),
    ("c1cc2cnncc2[nH]1",  "1H-pyrrolo[2,3-d]pyridazine"),
    ("c1cc2c[nH]cc2nn1",  "1H-pyrrolo[3,4-c]pyridazine"),
    ("c1cc2cc[nH]c2nn1",  "1H-pyrrolo[2,3-c]pyridazine"),
])
def test_phase494_pyridazine_nh_fusions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
