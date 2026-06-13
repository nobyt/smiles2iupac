"""Phase 477: isothiazolo[x,y-c]pyridine bicyclic retained names
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cc2csnc2cn1",  "isothiazolo[3,4-c]pyridine"),
    ("c1cc2nscc2cn1",  "isothiazolo[4,3-c]pyridine"),
    ("c1cc2sncc2cn1",  "isothiazolo[4,5-c]pyridine"),
    ("c1cc2cnsc2cn1",  "isothiazolo[5,4-c]pyridine"),
])
def test_phase477_isothiazolo_c_pyridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
