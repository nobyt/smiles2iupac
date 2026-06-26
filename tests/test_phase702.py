"""Phase 702: methyl derivatives of [1,2,4]triazolo[1,5-a] and [4,3-a/b] series
(9 parent compounds: pyridine, pyrimidine, pyrazine, pyridazine, [1,2,4]triazine).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # [1,2,4]triazolo[1,5-a]pyridine (CH at 2,5,6,7,8)
    ("c1ccn2ncnc2c1",      "[1,2,4]triazolo[1,5-a]pyridine"),
    ("Cc1nc2ccccn2n1",     "2-methyl[1,2,4]triazolo[1,5-a]pyridine"),
    ("Cc1cccc2ncnn12",     "5-methyl[1,2,4]triazolo[1,5-a]pyridine"),
    ("Cc1ccc2ncnn2c1",     "6-methyl[1,2,4]triazolo[1,5-a]pyridine"),
    ("Cc1ccn2ncnc2c1",     "7-methyl[1,2,4]triazolo[1,5-a]pyridine"),
    ("Cc1cccn2ncnc12",     "8-methyl[1,2,4]triazolo[1,5-a]pyridine"),
    # [1,2,4]triazolo[1,5-a]pyrimidine (CH at 2,5,6,7)
    ("c1cnc2ncnn2c1",      "[1,2,4]triazolo[1,5-a]pyrimidine"),
    ("Cc1nc2ncccn2n1",     "2-methyl[1,2,4]triazolo[1,5-a]pyrimidine"),
    ("Cc1ccn2ncnc2n1",     "5-methyl[1,2,4]triazolo[1,5-a]pyrimidine"),
    ("Cc1cnc2ncnn2c1",     "6-methyl[1,2,4]triazolo[1,5-a]pyrimidine"),
    ("Cc1ccnc2ncnn12",     "7-methyl[1,2,4]triazolo[1,5-a]pyrimidine"),
    # [1,2,4]triazolo[1,5-a]pyrazine (CH at 2,5,6,8)
    ("c1cn2ncnc2cn1",      "[1,2,4]triazolo[1,5-a]pyrazine"),
    ("Cc1nc2cnccn2n1",     "2-methyl[1,2,4]triazolo[1,5-a]pyrazine"),
    ("Cc1cncc2ncnn12",     "5-methyl[1,2,4]triazolo[1,5-a]pyrazine"),
    ("Cc1cn2ncnc2cn1",     "6-methyl[1,2,4]triazolo[1,5-a]pyrazine"),
    ("Cc1nccn2ncnc12",     "8-methyl[1,2,4]triazolo[1,5-a]pyrazine"),
    # [1,2,4]triazolo[1,5-b]pyridazine (CH at 2,6,7,8)
    ("c1cnn2ncnc2c1",      "[1,2,4]triazolo[1,5-b]pyridazine"),
    ("Cc1nc2cccnn2n1",     "2-methyl[1,2,4]triazolo[1,5-b]pyridazine"),
    ("Cc1ccc2ncnn2n1",     "6-methyl[1,2,4]triazolo[1,5-b]pyridazine"),
    ("Cc1cnn2ncnc2c1",     "7-methyl[1,2,4]triazolo[1,5-b]pyridazine"),
    ("Cc1ccnn2ncnc12",     "8-methyl[1,2,4]triazolo[1,5-b]pyridazine"),
    # [1,2,4]triazolo[1,5-b][1,2,4]triazine (CH at 2,6,7)
    ("c1cnn2ncnc2n1",      "[1,2,4]triazolo[1,5-b][1,2,4]triazine"),
    ("Cc1nc2nccnn2n1",     "2-methyl[1,2,4]triazolo[1,5-b][1,2,4]triazine"),
    ("Cc1cnc2ncnn2n1",     "6-methyl[1,2,4]triazolo[1,5-b][1,2,4]triazine"),
    ("Cc1cnn2ncnc2n1",     "7-methyl[1,2,4]triazolo[1,5-b][1,2,4]triazine"),
    # [1,2,4]triazolo[4,3-a]pyridine (CH at 3,5,6,7,8)
    ("c1ccn2cnnc2c1",      "[1,2,4]triazolo[4,3-a]pyridine"),
    ("Cc1nnc2ccccn12",     "3-methyl[1,2,4]triazolo[4,3-a]pyridine"),
    ("Cc1cccc2nncn12",     "5-methyl[1,2,4]triazolo[4,3-a]pyridine"),
    ("Cc1ccc2nncn2c1",     "6-methyl[1,2,4]triazolo[4,3-a]pyridine"),
    ("Cc1ccn2cnnc2c1",     "7-methyl[1,2,4]triazolo[4,3-a]pyridine"),
    ("Cc1cccn2cnnc12",     "8-methyl[1,2,4]triazolo[4,3-a]pyridine"),
    # [1,2,4]triazolo[4,3-a]pyrimidine (CH at 3,5,6,7)
    ("c1cnc2nncn2c1",      "[1,2,4]triazolo[4,3-a]pyrimidine"),
    ("Cc1nnc2ncccn12",     "3-methyl[1,2,4]triazolo[4,3-a]pyrimidine"),
    ("Cc1ccnc2nncn12",     "5-methyl[1,2,4]triazolo[4,3-a]pyrimidine"),
    ("Cc1cnc2nncn2c1",     "6-methyl[1,2,4]triazolo[4,3-a]pyrimidine"),
    ("Cc1ccn2cnnc2n1",     "7-methyl[1,2,4]triazolo[4,3-a]pyrimidine"),
    # [1,2,4]triazolo[4,3-a]pyrazine (CH at 3,5,6,8)
    ("c1cn2cnnc2cn1",      "[1,2,4]triazolo[4,3-a]pyrazine"),
    ("Cc1nnc2cnccn12",     "3-methyl[1,2,4]triazolo[4,3-a]pyrazine"),
    ("Cc1cncc2nncn12",     "5-methyl[1,2,4]triazolo[4,3-a]pyrazine"),
    ("Cc1cn2cnnc2cn1",     "6-methyl[1,2,4]triazolo[4,3-a]pyrazine"),
    ("Cc1nccn2cnnc12",     "8-methyl[1,2,4]triazolo[4,3-a]pyrazine"),
    # [1,2,4]triazolo[4,3-b][1,2,4]triazine (CH at 3,6,7)
    ("c1cnn2cnnc2n1",      "[1,2,4]triazolo[4,3-b][1,2,4]triazine"),
    ("Cc1nnc2nccnn12",     "3-methyl[1,2,4]triazolo[4,3-b][1,2,4]triazine"),
    ("Cc1cnc2nncn2n1",     "6-methyl[1,2,4]triazolo[4,3-b][1,2,4]triazine"),
    ("Cc1cnn2cnnc2n1",     "7-methyl[1,2,4]triazolo[4,3-b][1,2,4]triazine"),
])
def test_phase702(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
