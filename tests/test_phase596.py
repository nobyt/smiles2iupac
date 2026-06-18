"""Phase 596: Substituted imidazo[1,2-b]pyridazine, imidazo[1,5-b]pyridazine,
pyrazolo[1,5-b]pyridazine, [1,2,4]triazolo[1,5-b]pyridazine,
[1,2,3]triazolo[1,5-b]pyridazine, tetrazolo[1,5-b]pyridazine,
imidazo[1,2-a]pyrimidine, and imidazo[1,5-a]pyrimidine naming.
N-bridged 5+6 bicyclics with pyridazine or pyrimidine as the 6-ring.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imidazo[1,2-b]pyridazine (N at 4,9; sub C: 2,3,6,7,8)
    ("c1cnn2ccnc2c1",          "imidazo[1,2-b]pyridazine"),
    ("Cc1cn2ncccc2n1",         "2-methylimidazo[1,2-b]pyridazine"),
    ("Cc1cnc2cccnn12",         "3-methylimidazo[1,2-b]pyridazine"),
    ("Cc1ccc2nccn2n1",         "6-methylimidazo[1,2-b]pyridazine"),
    ("Cc1cnn2ccnc2c1",         "7-methylimidazo[1,2-b]pyridazine"),
    ("Cc1ccnn2ccnc12",         "8-methylimidazo[1,2-b]pyridazine"),
    # imidazo[1,5-b]pyridazine (N at 4,5,9; sub C: 2,3,5,7)
    ("c1cnn2cncc2c1",          "imidazo[1,5-b]pyridazine"),
    ("Cc1ccc2cncn2n1",         "2-methylimidazo[1,5-b]pyridazine"),
    ("Cc1cnn2cncc2c1",         "3-methylimidazo[1,5-b]pyridazine"),
    ("Cc1ncn2ncccc12",         "5-methylimidazo[1,5-b]pyridazine"),
    ("Cc1ncc2cccnn12",         "7-methylimidazo[1,5-b]pyridazine"),
    # pyrazolo[1,5-b]pyridazine (N at 4,5,9; sub C: 2,3,5,6)
    ("c1cnn2nccc2c1",          "pyrazolo[1,5-b]pyridazine"),
    ("Cc1cc2cccnn2n1",         "2-methylpyrazolo[1,5-b]pyridazine"),
    ("Cc1cnn2ncccc12",         "3-methylpyrazolo[1,5-b]pyridazine"),
    ("Cc1cnn2nccc2c1",         "5-methylpyrazolo[1,5-b]pyridazine"),
    ("Cc1ccc2ccnn2n1",         "6-methylpyrazolo[1,5-b]pyridazine"),
    # [1,2,4]triazolo[1,5-b]pyridazine (N at 4,5,6,9; sub C: 2,6,7,8)
    ("c1cnn2ncnc2c1",          "[1,2,4]triazolo[1,5-b]pyridazine"),
    ("Cc1nc2cccnn2n1",         "2-methyl[1,2,4]triazolo[1,5-b]pyridazine"),
    ("Cc1ccc2ncnn2n1",         "6-methyl[1,2,4]triazolo[1,5-b]pyridazine"),
    ("Cc1cnn2ncnc2c1",         "7-methyl[1,2,4]triazolo[1,5-b]pyridazine"),
    ("Cc1ccnn2ncnc12",         "8-methyl[1,2,4]triazolo[1,5-b]pyridazine"),
    # [1,2,3]triazolo[1,5-b]pyridazine (N at 4,5,6,9; sub C: 3,5,6)
    ("c1cnn2nncc2c1",          "[1,2,3]triazolo[1,5-b]pyridazine"),
    ("Cc1nnn2ncccc12",         "3-methyl[1,2,3]triazolo[1,5-b]pyridazine"),
    ("Cc1cnn2nncc2c1",         "5-methyl[1,2,3]triazolo[1,5-b]pyridazine"),
    ("Cc1ccc2cnnn2n1",         "6-methyl[1,2,3]triazolo[1,5-b]pyridazine"),
    # tetrazolo[1,5-b]pyridazine (N at 4,5,6,7,9; sub C: 6,7,8)
    ("c1cnn2nnnc2c1",          "tetrazolo[1,5-b]pyridazine"),
    ("Cc1ccc2nnnn2n1",         "6-methyltetrazolo[1,5-b]pyridazine"),
    ("Cc1cnn2nnnc2c1",         "7-methyltetrazolo[1,5-b]pyridazine"),
    ("Cc1ccnn2nnnc12",         "8-methyltetrazolo[1,5-b]pyridazine"),
    # imidazo[1,2-a]pyrimidine (N at 4,8,9; sub C: 2,3,5,6,7)
    ("c1cnc2nccn2c1",          "imidazo[1,2-a]pyrimidine"),
    ("Cc1cn2cccnc2n1",         "2-methylimidazo[1,2-a]pyrimidine"),
    ("Cc1cnc2ncccn12",         "3-methylimidazo[1,2-a]pyrimidine"),
    ("Cc1ccnc2nccn12",         "5-methylimidazo[1,2-a]pyrimidine"),
    ("Cc1cnc2nccn2c1",         "6-methylimidazo[1,2-a]pyrimidine"),
    ("Cc1ccn2ccnc2n1",         "7-methylimidazo[1,2-a]pyrimidine"),
    # imidazo[1,5-a]pyrimidine (N at 4,5,8,9; sub C: 2,3,6,8)
    ("c1cnc2cncn2c1",          "imidazo[1,5-a]pyrimidine"),
    ("Cc1ccn2cncc2n1",         "2-methylimidazo[1,5-a]pyrimidine"),
    ("Cc1cnc2cncn2c1",         "3-methylimidazo[1,5-a]pyrimidine"),
    ("Cc1ncc2ncccn12",         "6-methylimidazo[1,5-a]pyrimidine"),
    ("Cc1ncn2cccnc12",         "8-methylimidazo[1,5-a]pyrimidine"),
])
def test_phase596_n_bridged_pyridazine_pyrimidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
