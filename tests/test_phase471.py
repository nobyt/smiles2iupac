"""Phase 471: thieno/furo/pyrrolo[x,y-d]pyrimidine bicyclic retained names
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1ncc2ccsc2n1",    "thieno[2,3-d]pyrimidine"),
    ("c1ncc2cscc2n1",    "thieno[3,4-d]pyrimidine"),
    ("c1ncc2sccc2n1",    "thieno[3,2-d]pyrimidine"),
    ("c1ncc2ccoc2n1",    "furo[2,3-d]pyrimidine"),
    ("c1ncc2cocc2n1",    "furo[3,4-d]pyrimidine"),
    ("c1ncc2occc2n1",    "furo[3,2-d]pyrimidine"),
    ("c1ncc2cc[nH]c2n1", "1H-pyrrolo[2,3-d]pyrimidine"),
    ("c1ncc2c[nH]cc2n1", "1H-pyrrolo[3,4-d]pyrimidine"),
    ("c1ncc2[nH]ccc2n1", "1H-pyrrolo[3,2-d]pyrimidine"),
])
def test_phase471_thieno_furo_pyrrolo_pyrimidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
