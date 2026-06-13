"""Phase 478: isothiazolo[x,y-c]pyridazine bicyclic retained names
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2csnc2nn1",  "isothiazolo[3,4-c]pyridazine"),
    ("c1cc2nscc2nn1",  "isothiazolo[4,3-c]pyridazine"),
    ("c1cc2sncc2nn1",  "isothiazolo[4,5-c]pyridazine"),
    ("c1cc2cnsc2nn1",  "isothiazolo[5,4-c]pyridazine"),
])
def test_phase478_isothiazolo_c_pyridazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
