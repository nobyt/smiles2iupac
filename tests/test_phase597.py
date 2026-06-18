"""Phase 597: Substituted pyrazolo[1,5-a]pyrimidine, [1,2,4]triazolo[4,3-a]pyrimidine,
[1,2,4]triazolo[1,5-a]pyrimidine, pyrrolo[1,2-a]pyrazine,
imidazo[1,2-a]pyrazine, pyrrolo[1,2-b]pyridazine, and
pyrrolo[1,2-a]pyrimidine naming.
N-bridged 5+6 bicyclics with pyrimidine or pyrazine as the 6-ring.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrazolo[1,5-a]pyrimidine (N at 4,5; sub C: 2,3,5,6,7)
    ("c1cnc2ccnn2c1",          "pyrazolo[1,5-a]pyrimidine"),
    ("Cc1cc2ncccn2n1",         "2-methylpyrazolo[1,5-a]pyrimidine"),
    ("Cc1cnn2cccnc12",         "3-methylpyrazolo[1,5-a]pyrimidine"),
    ("Cc1ccn2nccc2n1",         "5-methylpyrazolo[1,5-a]pyrimidine"),
    ("Cc1cnc2ccnn2c1",         "6-methylpyrazolo[1,5-a]pyrimidine"),
    ("Cc1ccnc2ccnn12",         "7-methylpyrazolo[1,5-a]pyrimidine"),
    # [1,2,4]triazolo[4,3-a]pyrimidine (N at 4,6,7; sub C: 3,5,6,7)
    ("c1cnc2nncn2c1",          "[1,2,4]triazolo[4,3-a]pyrimidine"),
    ("Cc1nnc2ncccn12",         "3-methyl[1,2,4]triazolo[4,3-a]pyrimidine"),
    ("Cc1ccnc2nncn12",         "5-methyl[1,2,4]triazolo[4,3-a]pyrimidine"),
    ("Cc1cnc2nncn2c1",         "6-methyl[1,2,4]triazolo[4,3-a]pyrimidine"),
    ("Cc1ccn2cnnc2n1",         "7-methyl[1,2,4]triazolo[4,3-a]pyrimidine"),
    # [1,2,4]triazolo[1,5-a]pyrimidine (N at 4,5,7; sub C: 2,5,6,7)
    ("c1cnc2ncnn2c1",          "[1,2,4]triazolo[1,5-a]pyrimidine"),
    ("Cc1nc2ncccn2n1",         "2-methyl[1,2,4]triazolo[1,5-a]pyrimidine"),
    ("Cc1ccn2ncnc2n1",         "5-methyl[1,2,4]triazolo[1,5-a]pyrimidine"),
    ("Cc1cnc2ncnn2c1",         "6-methyl[1,2,4]triazolo[1,5-a]pyrimidine"),
    ("Cc1ccnc2ncnn12",         "7-methyl[1,2,4]triazolo[1,5-a]pyrimidine"),
    # pyrrolo[1,2-a]pyrazine (N at 4,9; sub C: 1,3,4,6,7,8)
    ("c1cc2cnccn2c1",          "pyrrolo[1,2-a]pyrazine"),
    ("Cc1nccn2cccc12",         "1-methylpyrrolo[1,2-a]pyrazine"),
    ("Cc1cn2cccc2cn1",         "3-methylpyrrolo[1,2-a]pyrazine"),
    ("Cc1cncc2cccn12",         "4-methylpyrrolo[1,2-a]pyrazine"),
    ("Cc1ccc2cnccn12",         "6-methylpyrrolo[1,2-a]pyrazine"),
    ("Cc1cc2cnccn2c1",         "7-methylpyrrolo[1,2-a]pyrazine"),
    ("Cc1ccn2ccncc12",         "8-methylpyrrolo[1,2-a]pyrazine"),
    # imidazo[1,2-a]pyrazine (N at 4,5,9; sub C: 2,3,5,6,8)
    ("c1cn2ccnc2cn1",          "imidazo[1,2-a]pyrazine"),
    ("Cc1cn2ccncc2n1",         "2-methylimidazo[1,2-a]pyrazine"),
    ("Cc1cnc2cnccn12",         "3-methylimidazo[1,2-a]pyrazine"),
    ("Cc1cncc2nccn12",         "5-methylimidazo[1,2-a]pyrazine"),
    ("Cc1cn2ccnc2cn1",         "6-methylimidazo[1,2-a]pyrazine"),
    ("Cc1nccn2ccnc12",         "8-methylimidazo[1,2-a]pyrazine"),
    # pyrrolo[1,2-b]pyridazine (N at 1,2,9; sub C: 2,3,4,5,6,7)
    ("c1cnn2cccc2c1",          "pyrrolo[1,2-b]pyridazine"),
    ("Cc1ccc2cccn2n1",         "2-methylpyrrolo[1,2-b]pyridazine"),
    ("Cc1cnn2cccc2c1",         "3-methylpyrrolo[1,2-b]pyridazine"),
    ("Cc1ccnn2cccc12",         "4-methylpyrrolo[1,2-b]pyridazine"),
    ("Cc1ccn2ncccc12",         "5-methylpyrrolo[1,2-b]pyridazine"),
    ("Cc1cc2cccnn2c1",         "6-methylpyrrolo[1,2-b]pyridazine"),
    ("Cc1ccc2cccnn12",         "7-methylpyrrolo[1,2-b]pyridazine"),
    # pyrrolo[1,2-a]pyrimidine (N at 4,9; sub C: 2,3,4,6,7,8)
    ("c1cnc2cccn2c1",          "pyrrolo[1,2-a]pyrimidine"),
    ("Cc1ccn2cccc2n1",         "2-methylpyrrolo[1,2-a]pyrimidine"),
    ("Cc1cnc2cccn2c1",         "3-methylpyrrolo[1,2-a]pyrimidine"),
    ("Cc1ccnc2cccn12",         "4-methylpyrrolo[1,2-a]pyrimidine"),
    ("Cc1ccc2ncccn12",         "6-methylpyrrolo[1,2-a]pyrimidine"),
    ("Cc1cc2ncccn2c1",         "7-methylpyrrolo[1,2-a]pyrimidine"),
    ("Cc1ccn2cccnc12",         "8-methylpyrrolo[1,2-a]pyrimidine"),
])
def test_phase597_n_bridged_pyrimidine_pyrazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
