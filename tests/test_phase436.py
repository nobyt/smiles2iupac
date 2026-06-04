"""Phase 436: More drug-scaffold bicyclic heterocycles — furo-pyrimidine,
pyrazolo-pyridine, triazolo-pyridine (IUPAC 2013 P-31.1.3).

These bicyclics currently output 'benzene' (wrong).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # furo[2,3-d]pyrimidine — O adjacent to ring junction (kinase scaffold)
    ("c1cnc2occc2n1",              "furo[2,3-d]pyrimidine"),
    # furo[3,2-d]pyrimidine — O at mid-position of furan ring
    ("c1cnc2cocc2n1",              "furo[3,2-d]pyrimidine"),
    # pyrazolo[1,5-a]pyridine
    ("c1ccn2nccc2c1",              "pyrazolo[1,5-a]pyridine"),
    # [1,2,4]triazolo[4,3-a]pyridine
    ("c1ccn2cnnc2c1",              "[1,2,4]triazolo[4,3-a]pyridine"),
    # furo[2,3-b]pyridine
    ("c1cnc2occc2c1",              "furo[2,3-b]pyridine"),
    # 1H-pyrazolo[3,4-b]pyridine
    ("c1cnc2[nH]ncc2c1",           "1H-pyrazolo[3,4-b]pyridine"),
    # regression: imidazo[1,2-a]pyridine unchanged
    ("c1ccn2ccnc2c1",              "imidazo[1,2-a]pyridine"),
    # regression: 7H-pyrrolo[2,3-d]pyrimidine unchanged (Phase 435)
    ("c1cnc2[nH]ccc2n1",           "7H-pyrrolo[2,3-d]pyrimidine"),
])
def test_phase436_more_drug_scaffolds(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
