"""Phase 442: imidazo[1,5-a]pyridine, tetrazolo-pyrimidine, pyrrolo-pyridazine,
thieno/furo-pyridine isomers, and 1H-pyrazolo[3,4-d]pyrimidine retained names
(IUPAC 2013 P-31.1.3).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imidazo[1,5-a]pyridine
    ("c1ccn2cncc2c1",       "imidazo[1,5-a]pyridine"),
    # tetrazolo[1,5-a]pyrimidine
    ("c1cnc2nnnn2c1",       "tetrazolo[1,5-a]pyrimidine"),
    # pyrrolo[1,2-a]pyrimidine
    ("c1cnc2cccn2c1",       "pyrrolo[1,2-a]pyrimidine"),
    # thieno[3,2-c]pyridine
    ("c1cc2sccc2cn1",       "thieno[3,2-c]pyridine"),
    # furo[3,2-b]pyridine
    ("c1cnc2ccoc2c1",       "furo[3,2-b]pyridine"),
    # furo[3,4-b]pyridine
    ("c1cnc2cocc2c1",       "furo[3,4-b]pyridine"),
    # 1H-pyrazolo[3,4-d]pyrimidine
    ("c1ncc2c[nH]nc2n1",  "1H-pyrazolo[3,4-d]pyrimidine"),
    # regressions
    ("c1ccn2cncc2c1",       "imidazo[1,5-a]pyridine"),
    ("c1cnc2occc2c1",       "furo[2,3-b]pyridine"),
    ("c1cnc2ncnn2c1",       "[1,2,4]triazolo[1,5-a]pyrimidine"),
    ("c1cc2nccnn2c1",       "pyrrolo[1,2-b][1,2,4]triazine"),
    ("c1ccn2nnnc2c1",       "tetrazolo[1,5-a]pyridine"),
])
def test_phase442_fused_hetero(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
