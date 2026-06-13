"""Phase 476: isothiazolo[x,y-b]pyridine bicyclic retained names
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 4 C-C bond junctions of isothiazole with pyridine b-bond (C2-C3)
    ("c1cnc2nscc2c1",  "isothiazolo[3,4-b]pyridine"),
    ("c1cnc2csnc2c1",  "isothiazolo[4,3-b]pyridine"),
    ("c1cnc2cnsc2c1",  "isothiazolo[4,5-b]pyridine"),
    ("c1cnc2sncc2c1",  "isothiazolo[5,4-b]pyridine"),
])
def test_phase476_isothiazolo_b_pyridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
