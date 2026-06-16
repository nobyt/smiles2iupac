"""Phase 435: Drug-scaffold bicyclic heterocycles (IUPAC 2013 P-31.1.3).

Names corrected in Phase 473: c1cnc2 ring6 is pyrazine (N at 1,4), not pyrimidine.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-pyrrolo[2,3-e]pyrazine (corrected Phase 473: ring6 is pyrazine)
    ("c1c[nH]c2ccnc-2n1",           "1H-pyrrolo[2,3-e]pyrazine"),
    # 1H-pyrazolo[4,5-e]pyrazine (corrected Phase 473/474: [4,5]<[5,4] preferred)
    ("c1cnc2[nH]ncc2n1",           "1H-pyrazolo[4,5-e]pyrazine"),
    # 1H-imidazo[4,5-c]pyridine
    ("c1cc2nc[nH]c2cn1",           "1H-imidazo[4,5-c]pyridine"),
    # thieno[2,3-e]pyrazine (corrected Phase 473)
    ("c1cnc2sccc2n1",              "thieno[2,3-e]pyrazine"),
    # regression: 1H-pyrrolo[2,3-b]pyridine unchanged (Phase 255)
    ("c1cnc2[nH]ccc2c1",           "1H-pyrrolo[2,3-b]pyridine"),
    # regression: 9H-purine (NH at N9, adjacent to C4)
    ("c1ncc2nc[nH]c2n1",           "9H-purine"),
    # regression: pyrimidine unchanged
    ("c1ccncn1",                   "pyrimidine"),
])
def test_phase435_drug_scaffold_bicyclics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
