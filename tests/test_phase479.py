"""Phase 479: missing isoxazolo isomers completing the [x,y-b/c] families
(IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isoxazolo[4,3-b]pyridine (jA=iso-C4, jB=iso-C3 adj-N, completing [3,4]/[4,5]/[5,4] set)
    ("c1cnc2conc2c1",  "isoxazolo[4,3-b]pyridine"),
    # isoxazolo[3,4-c]pyridine (jA=iso-C4, jB=iso-C3 adj-N, completing [4,3]/[4,5]/[5,4] set)
    ("c1cc2conc2cn1",  "isoxazolo[3,4-c]pyridine"),
    # isoxazolo[3,4-c]pyridazine (completing the c-pyridazine set)
    ("c1cc2conc2nn1",  "isoxazolo[3,4-c]pyridazine"),
])
def test_phase479_isoxazolo_missing(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
