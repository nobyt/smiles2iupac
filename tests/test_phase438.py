"""Phase 438: pyrido-pyrimidine, pyrimido-pyrimidine, and thieno-pyridine
bicyclic retained names (IUPAC 2013 P-31.1.3).

These bicyclics currently output 'benzene' (wrong).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrido[2,3-d]pyrimidine — C7H5N3 bicyclic
    ("c1cnc2nccnc2c1",             "pyrido[2,3-d]pyrimidine"),
    # pyrido[3,4-d]pyrimidine — different N arrangement
    ("c1cnc2cncnc2c1",             "pyrido[3,4-d]pyrimidine"),
    # pyrimido[4,5-d]pyrimidine — C6H4N4 bicyclic (all-N-substituted positions)
    ("c1cnc2ncncc2n1",             "pyrimido[4,5-d]pyrimidine"),
    # thieno[3,4-b]pyridine — C7H5NS bicyclic
    ("c1cnc2cscc2c1",              "thieno[3,4-b]pyridine"),
    # regression: pteridine unchanged (correct output already)
    ("c1cnc2nccnc2n1",             "pteridine"),
    # regression: 6H-purine unchanged
    ("c1ncc2nc[nH]c2n1",           "6H-purine"),
])
def test_phase438_pyrido_pyrimidines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
