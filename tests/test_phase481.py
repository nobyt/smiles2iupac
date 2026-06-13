"""Phase 481: imidazo[3,2-b]pyridazine missing fused name
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cnn2ccnc2c1",  "imidazo[1,2-b]pyridazine"),
])
def test_phase481_imidazo_pyridazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
