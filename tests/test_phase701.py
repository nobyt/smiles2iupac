"""Phase 701: methyl derivatives of imidazo[3,2-b]/[1,5-b][1,2,4]triazine
and pyrazolo[1,5-a/b] series (7 parent compounds).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imidazo[3,2-b][1,2,4]triazine (CH at 2,3,6,7)
    ("c1cnn2ccnc2n1",      "imidazo[3,2-b][1,2,4]triazine"),
    ("Cc1cnc2nccn2n1",     "2-methylimidazo[3,2-b][1,2,4]triazine"),
    ("Cc1cnn2ccnc2n1",     "3-methylimidazo[3,2-b][1,2,4]triazine"),
    ("Cc1cn2nccnc2n1",     "6-methylimidazo[3,2-b][1,2,4]triazine"),
    ("Cc1cnc2nccnn12",     "7-methylimidazo[3,2-b][1,2,4]triazine"),
    # imidazo[1,5-b][1,2,4]triazine (CH at 2,3,6,8)
    ("c1cnn2cncc2n1",      "imidazo[1,5-b][1,2,4]triazine"),
    ("Cc1cnn2cncc2n1",     "2-methylimidazo[1,5-b][1,2,4]triazine"),
    ("Cc1cnc2cncn2n1",     "3-methylimidazo[1,5-b][1,2,4]triazine"),
    ("Cc1ncc2nccnn12",     "6-methylimidazo[1,5-b][1,2,4]triazine"),
    ("Cc1ncn2nccnc12",     "8-methylimidazo[1,5-b][1,2,4]triazine"),
    # pyrazolo[1,5-a]pyridine (CH at 2,3,5,6,7)
    ("c1ccn2nccc2c1",      "pyrazolo[1,5-a]pyridine"),
    ("Cc1cc2ccccn2n1",     "2-methylpyrazolo[1,5-a]pyridine"),
    ("Cc1cnn2ccccc12",     "3-methylpyrazolo[1,5-a]pyridine"),
    ("Cc1ccn2nccc2c1",     "5-methylpyrazolo[1,5-a]pyridine"),
    ("Cc1ccc2ccnn2c1",     "6-methylpyrazolo[1,5-a]pyridine"),
    ("Cc1cccc2ccnn12",     "7-methylpyrazolo[1,5-a]pyridine"),
    # pyrazolo[1,5-a]pyrimidine (CH at 2,3,5,6,7)
    ("c1cnc2ccnn2c1",      "pyrazolo[1,5-a]pyrimidine"),
    ("Cc1cc2ncccn2n1",     "2-methylpyrazolo[1,5-a]pyrimidine"),
    ("Cc1cnn2cccnc12",     "3-methylpyrazolo[1,5-a]pyrimidine"),
    ("Cc1ccn2nccc2n1",     "5-methylpyrazolo[1,5-a]pyrimidine"),
    ("Cc1cnc2ccnn2c1",     "6-methylpyrazolo[1,5-a]pyrimidine"),
    ("Cc1ccnc2ccnn12",     "7-methylpyrazolo[1,5-a]pyrimidine"),
    # pyrazolo[1,5-a]pyrazine (CH at 2,3,4,6,7)
    ("c1cn2nccc2cn1",      "pyrazolo[1,5-a]pyrazine"),
    ("Cc1cc2cnccn2n1",     "2-methylpyrazolo[1,5-a]pyrazine"),
    ("Cc1cnn2ccncc12",     "3-methylpyrazolo[1,5-a]pyrazine"),
    ("Cc1nccn2nccc12",     "4-methylpyrazolo[1,5-a]pyrazine"),
    ("Cc1cn2nccc2cn1",     "6-methylpyrazolo[1,5-a]pyrazine"),
    ("Cc1cncc2ccnn12",     "7-methylpyrazolo[1,5-a]pyrazine"),
    # pyrazolo[1,5-b]pyridazine (CH at 2,3,5,6)
    ("c1cnn2nccc2c1",      "pyrazolo[1,5-b]pyridazine"),
    ("Cc1cc2cccnn2n1",     "2-methylpyrazolo[1,5-b]pyridazine"),
    ("Cc1cnn2ncccc12",     "3-methylpyrazolo[1,5-b]pyridazine"),
    ("Cc1cnn2nccc2c1",     "5-methylpyrazolo[1,5-b]pyridazine"),
    ("Cc1ccc2ccnn2n1",     "6-methylpyrazolo[1,5-b]pyridazine"),
    # pyrazolo[1,5-b][1,2,4]triazine (CH at 2,3,7,8)
    ("c1cnn2nccc2n1",      "pyrazolo[1,5-b][1,2,4]triazine"),
    ("Cc1cnn2nccc2n1",     "2-methylpyrazolo[1,5-b][1,2,4]triazine"),
    ("Cc1cnc2ccnn2n1",     "3-methylpyrazolo[1,5-b][1,2,4]triazine"),
    ("Cc1cc2nccnn2n1",     "7-methylpyrazolo[1,5-b][1,2,4]triazine"),
    ("Cc1cnn2nccnc12",     "8-methylpyrazolo[1,5-b][1,2,4]triazine"),
])
def test_phase701(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
