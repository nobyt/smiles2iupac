"""Phase 436: More drug-scaffold bicyclic heterocycles (IUPAC 2013 P-31.1.3).

furo-pyrimidine entries corrected in Phase 473 (ring6 is pyrazine, not pyrimidine).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # furo[2,3-e]pyrazine (corrected Phase 473: ring6 is pyrazine)
    ("c1cnc2occc2n1",              "furo[2,3-e]pyrazine"),
    # furo[3,4-e]pyrazine (corrected Phase 473)
    ("c1cnc2cocc2n1",              "furo[3,4-e]pyrazine"),
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
    # regression: 1H-pyrrolo[2,3-e]pyrazine unchanged (Phase 435, corrected Phase 473)
    ("c1cnc2[nH]ccc2n1",           "1H-pyrrolo[2,3-e]pyrazine"),
])
def test_phase436_more_drug_scaffolds(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
