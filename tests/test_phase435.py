"""Phase 435: Drug-scaffold bicyclic heterocycles — pyrrolo/pyrazolo-pyrimidine and
imidazo/thieno-pyridine systems (IUPAC 2013 P-31.1.3).

These C6H5N3/C5H4N4/C6H4N2S bicyclics currently output 'benzene' (wrong).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 7H-pyrrolo[2,3-d]pyrimidine (beta-carboline analog, common kinase scaffold)
    ("c1cnc2[nH]ccc2n1",           "7H-pyrrolo[2,3-d]pyrimidine"),
    # 1H-pyrazolo[3,4-d]pyrimidine (adenine isostere scaffold)
    ("c1cnc2[nH]ncc2n1",           "1H-pyrazolo[3,4-d]pyrimidine"),
    # 1H-imidazo[4,5-c]pyridine
    ("c1cc2nc[nH]c2cn1",           "1H-imidazo[4,5-c]pyridine"),
    # thieno[3,2-d]pyrimidine (kinase inhibitor scaffold)
    ("c1cnc2sccc2n1",              "thieno[3,2-d]pyrimidine"),
    # regression: 1H-pyrrolo[2,3-b]pyridine unchanged (Phase 255)
    ("c1cnc2[nH]ccc2c1",           "1H-pyrrolo[2,3-b]pyridine"),
    # regression: 6H-purine unchanged
    ("c1ncc2nc[nH]c2n1",           "6H-purine"),
    # regression: pyrimidine unchanged
    ("c1ccncn1",                   "pyrimidine"),
])
def test_phase435_drug_scaffold_bicyclics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
