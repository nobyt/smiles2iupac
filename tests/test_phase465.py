"""Phase 465: furo[3,2-c]pyridine, furo[2,3-e]pyridazine,
1H-pyrazolo[3,4-c]pyridine, 1H-pyrrolo[3,2-c]pyridine
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2occc2cn1",    "furo[3,2-c]pyridine"),
    ("c1cc2occc2nn1",    "furo[2,3-e]pyridazine"),
    ("c1cc2[nH]ncc2cn1", "1H-pyrazolo[4,5-c]pyridine"),
    ("c1cc2[nH]ccc2cn1", "1H-pyrrolo[3,2-c]pyridine"),
])
def test_phase465_furo_pyrrolo_pyridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
