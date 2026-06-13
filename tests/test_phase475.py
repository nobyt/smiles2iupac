"""Phase 475: isothiazolo[x,y-d]pyrimidine and [x,y-e]pyrazine bicyclic retained names
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isothiazolo-d-pyrimidine (4 isomers: C3-C4 and C4-C5 bonds of isothiazole)
    ("c1ncc2sncc2n1",  "isothiazolo[4,5-d]pyrimidine"),
    ("c1ncc2nscc2n1",  "isothiazolo[4,3-d]pyrimidine"),
    ("c1ncc2cnsc2n1",  "isothiazolo[5,4-d]pyrimidine"),
    ("c1ncc2csnc2n1",  "isothiazolo[3,4-d]pyrimidine"),
    # isothiazolo-e-pyrazine (2 unique due to pyrazine C2 symmetry)
    ("c1cnc2sncc2n1",  "isothiazolo[4,5-e]pyrazine"),
    ("c1cnc2nscc2n1",  "isothiazolo[3,4-e]pyrazine"),
])
def test_phase475_isothiazolo_pyrimidine_pyrazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
