"""Phase 441: triazolo-pyrimidine, triazolo-pyridazine, pyrazolo-pyrimidine,
tetrazolo-pyridine, and 1H-imidazo[4,5-b]pyridine retained names
(IUPAC 2013 P-31.1.3).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # [1,2,4]triazolo[1,5-a]pyrimidine
    ("c1cnc2ncnn2c1",     "[1,2,4]triazolo[1,5-a]pyrimidine"),
    # [1,2,4]triazolo[1,5-b][1,2,4]triazine
    ("c1cnn2ncnc2n1",     "[1,2,4]triazolo[1,5-b][1,2,4]triazine"),
    # [1,2,4]triazolo[4,3-a]pyrimidine
    ("c1cnc2nncn2c1",     "[1,2,4]triazolo[4,3-a]pyrimidine"),
    # pyrrolo[1,2-b][1,2,4]triazine
    ("c1cc2nccnn2c1",     "pyrrolo[1,2-b][1,2,4]triazine"),
    # tetrazolo[1,5-a]pyridine
    ("c1ccn2nnnc2c1",     "tetrazolo[1,5-a]pyridine"),
    # 1H-imidazo[4,5-b]pyridine
    ("c1cnc2[nH]cnc2c1",  "1H-imidazo[4,5-b]pyridine"),
    # regressions
    ("c1ccn2cnnc2c1",     "[1,2,4]triazolo[4,3-a]pyridine"),
    ("c1cnc2cnnn2c1",     "[1,2,3]triazolo[1,5-a]pyrimidine"),
    ("c1cnc2[nH]nnc2c1",  "1H-[1,2,3]triazolo[4,5-b]pyridine"),
    ("c1cnc2[nH]cnc2n1",  "1H-imidazo[4,5-e]pyrazine"),  # corrected Phase 473
    ("c1cnc2nc[nH]c2c1",  "3H-imidazo[4,5-b]pyridine"),
])
def test_phase441_triazolo_pyrazolo_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
