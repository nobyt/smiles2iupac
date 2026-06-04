"""Phase 440: thieno-pyridine and triazolo-pyridine bicyclic retained names
(IUPAC 2013 P-31.1.3).

These bicyclics currently output 'benzene' (wrong).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[3,4-c]pyridine — S adjacent to ring junction
    ("c1cnc2sccc2c1",              "thieno[3,4-c]pyridine"),
    # [1,2,3]triazolo[1,5-a]pyridine
    ("c1ccn2nncc2c1",              "[1,2,3]triazolo[1,5-a]pyridine"),
    # 1H-[1,2,3]triazolo[4,5-b]pyridine
    ("c1cnc2[nH]nnc2c1",           "1H-[1,2,3]triazolo[4,5-b]pyridine"),
    # regression: thieno[3,4-b]pyridine unchanged (Phase 438)
    ("c1cnc2cscc2c1",              "thieno[3,4-b]pyridine"),
    # regression: [1,2,4]triazolo[4,3-a]pyridine unchanged (Phase 436)
    ("c1ccn2cnnc2c1",              "[1,2,4]triazolo[4,3-a]pyridine"),
    # regression: [1,2,3]triazolo[1,5-a]pyrimidine unchanged (Phase 437)
    ("c1cnc2cnnn2c1",              "[1,2,3]triazolo[1,5-a]pyrimidine"),
])
def test_phase440_thieno_triazolo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
