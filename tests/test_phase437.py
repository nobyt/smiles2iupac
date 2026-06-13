"""Phase 437: imidazo/thieno-pyrimidine and triazolo-pyrimidine bicyclics
(IUPAC 2013 P-31.1.3).

imidazo/thieno entries corrected in Phase 473 (c1cnc2 ring6 is pyrazine, not pyrimidine).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imidazo[1,2-a]pyrimidine
    ("c1cnc2nccn2c1",              "imidazo[1,2-a]pyrimidine"),
    # 1H-imidazo[4,5-e]pyrazine (corrected Phase 473: ring6 is pyrazine; [4,5]<[5,4] preferred)
    ("c1cnc2[nH]cnc2n1",           "1H-imidazo[4,5-e]pyrazine"),
    # thieno[3,4-e]pyrazine (corrected Phase 473)
    ("c1cnc2cscc2n1",              "thieno[3,4-e]pyrazine"),
    # [1,2,3]triazolo[1,5-a]pyrimidine
    ("c1cnc2cnnn2c1",              "[1,2,3]triazolo[1,5-a]pyrimidine"),
    # regression: imidazo[1,2-a]pyridine unchanged (Phase 255)
    ("c1ccn2ccnc2c1",              "imidazo[1,2-a]pyridine"),
    # regression: thieno[2,3-e]pyrazine unchanged (Phase 435, corrected Phase 473)
    ("c1cnc2sccc2n1",              "thieno[2,3-e]pyrazine"),
    # regression: [1,2,4]triazolo[4,3-a]pyridine unchanged (Phase 436)
    ("c1ccn2cnnc2c1",              "[1,2,4]triazolo[4,3-a]pyridine"),
])
def test_phase437_more_bicyclics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
