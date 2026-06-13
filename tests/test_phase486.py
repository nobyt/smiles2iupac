"""Phase 486: imidazo[4,5-c]pyridazine, thieno/furo-pyridazine isomers,
and pyrazolo[3,4-c]pyridine (IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2nc[nH]c2nn1",  "3H-imidazo[4,5-c]pyridazine"),
    ("c1cc2ccsc2nn1",     "thieno[2,3-c]pyridazine"),
    ("c1cc2ccoc2nn1",     "furo[2,3-c]pyridazine"),
    ("c1cc2n[nH]cc2cn1",  "1H-pyrazolo[3,4-c]pyridine"),
])
def test_phase486_fused_pyridazine_pyridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
