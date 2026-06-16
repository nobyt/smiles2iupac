"""Phase 438: pyrido-pyrimidine, pyrimido-pyrimidine, and thieno-pyridine
bicyclic retained names (IUPAC 2013 P-31.1.3).

These bicyclics currently output 'benzene' (wrong).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrido[2,3-d]pyrimidine — C7H5N3 bicyclic
    ("c1cnc2ncncc2c1",             "pyrido[2,3-d]pyrimidine"),
    # pyrido[3,4-d]pyrimidine — different N arrangement
    ("c1cc2cncnc2cn1",             "pyrido[3,4-d]pyrimidine"),
    # pteridine = pyrimido[4,5-d]pyrimidine — IUPAC retained name (P-31.1.3.4)
    ("c1cnc2ncncc2n1",             "pteridine"),
    # thieno[3,4-b]pyridine — C7H5NS bicyclic
    ("c1cnc2cscc2c1",              "thieno[3,4-b]pyridine"),
    # regression: 9H-purine (NH at N9, adjacent to C4)
    ("c1ncc2nc[nH]c2n1",           "9H-purine"),
])
def test_phase438_pyrido_pyrimidines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
