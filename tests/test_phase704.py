"""Phase 704: methyl derivatives of indolizine and pyrrolo[1,2-a/b] series
(5 parent compounds).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # indolizine (CH at 1,2,3,5,6,7,8)
    ("c1ccn2cccc2c1",      "indolizine"),
    ("Cc1ccn2ccccc12",     "1-methylindolizine"),
    ("Cc1cc2ccccn2c1",     "2-methylindolizine"),
    ("Cc1ccc2ccccn12",     "3-methylindolizine"),
    ("Cc1cccc2cccn12",     "5-methylindolizine"),
    ("Cc1ccc2cccn2c1",     "6-methylindolizine"),
    ("Cc1ccn2cccc2c1",     "7-methylindolizine"),
    ("Cc1cccn2cccc12",     "8-methylindolizine"),
    # pyrrolo[1,2-a]pyrazine (CH at 1,3,4,6,7,8)
    ("c1cc2cnccn2c1",      "pyrrolo[1,2-a]pyrazine"),
    ("Cc1nccn2cccc12",     "1-methylpyrrolo[1,2-a]pyrazine"),
    ("Cc1cn2cccc2cn1",     "3-methylpyrrolo[1,2-a]pyrazine"),
    ("Cc1cncc2cccn12",     "4-methylpyrrolo[1,2-a]pyrazine"),
    ("Cc1ccc2cnccn12",     "6-methylpyrrolo[1,2-a]pyrazine"),
    ("Cc1cc2cnccn2c1",     "7-methylpyrrolo[1,2-a]pyrazine"),
    ("Cc1ccn2ccncc12",     "8-methylpyrrolo[1,2-a]pyrazine"),
    # pyrrolo[1,2-a]pyrimidine (CH at 2,3,4,6,7,8)
    ("c1cnc2cccn2c1",      "pyrrolo[1,2-a]pyrimidine"),
    ("Cc1ccn2cccc2n1",     "2-methylpyrrolo[1,2-a]pyrimidine"),
    ("Cc1cnc2cccn2c1",     "3-methylpyrrolo[1,2-a]pyrimidine"),
    ("Cc1ccnc2cccn12",     "4-methylpyrrolo[1,2-a]pyrimidine"),
    ("Cc1ccc2ncccn12",     "6-methylpyrrolo[1,2-a]pyrimidine"),
    ("Cc1cc2ncccn2c1",     "7-methylpyrrolo[1,2-a]pyrimidine"),
    ("Cc1ccn2cccnc12",     "8-methylpyrrolo[1,2-a]pyrimidine"),
    # pyrrolo[1,2-b]pyridazine (CH at 2,3,4,5,6,7)
    ("c1cnn2cccc2c1",      "pyrrolo[1,2-b]pyridazine"),
    ("Cc1ccc2cccn2n1",     "2-methylpyrrolo[1,2-b]pyridazine"),
    ("Cc1cnn2cccc2c1",     "3-methylpyrrolo[1,2-b]pyridazine"),
    ("Cc1ccnn2cccc12",     "4-methylpyrrolo[1,2-b]pyridazine"),
    ("Cc1ccn2ncccc12",     "5-methylpyrrolo[1,2-b]pyridazine"),
    ("Cc1cc2cccnn2c1",     "6-methylpyrrolo[1,2-b]pyridazine"),
    ("Cc1ccc2cccnn12",     "7-methylpyrrolo[1,2-b]pyridazine"),
    # pyrrolo[1,2-b][1,2,4]triazine (CH at 2,3,6,7,8)
    ("c1cc2nccnn2c1",      "pyrrolo[1,2-b][1,2,4]triazine"),
    ("Cc1cnn2cccc2n1",     "2-methylpyrrolo[1,2-b][1,2,4]triazine"),
    ("Cc1cnc2cccn2n1",     "3-methylpyrrolo[1,2-b][1,2,4]triazine"),
    ("Cc1ccc2nccnn12",     "6-methylpyrrolo[1,2-b][1,2,4]triazine"),
    ("Cc1cc2nccnn2c1",     "7-methylpyrrolo[1,2-b][1,2,4]triazine"),
    ("Cc1ccn2nccnc12",     "8-methylpyrrolo[1,2-b][1,2,4]triazine"),
])
def test_phase704(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
