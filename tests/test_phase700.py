"""Phase 700: methyl derivatives of imidazo[1,2-a/b/c] and imidazo[1,5-a/b]
bicyclic series (9 parent compounds).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imidazo[1,2-a]pyridine (CH at 2,3,5,6,7,8)
    ("c1ccn2ccnc2c1",       "imidazo[1,2-a]pyridine"),
    ("Cc1cn2ccccc2n1",      "2-methylimidazo[1,2-a]pyridine"),
    ("Cc1cnc2ccccn12",      "3-methylimidazo[1,2-a]pyridine"),
    ("Cc1cccc2nccn12",      "5-methylimidazo[1,2-a]pyridine"),
    ("Cc1ccc2nccn2c1",      "6-methylimidazo[1,2-a]pyridine"),
    ("Cc1ccn2ccnc2c1",      "7-methylimidazo[1,2-a]pyridine"),
    ("Cc1cccn2ccnc12",      "8-methylimidazo[1,2-a]pyridine"),
    # imidazo[1,2-a]pyrimidine (CH at 2,3,5,6,7)
    ("c1cnc2nccn2c1",       "imidazo[1,2-a]pyrimidine"),
    ("Cc1cn2cccnc2n1",      "2-methylimidazo[1,2-a]pyrimidine"),
    ("Cc1cnc2ncccn12",      "3-methylimidazo[1,2-a]pyrimidine"),
    ("Cc1ccnc2nccn12",      "5-methylimidazo[1,2-a]pyrimidine"),
    ("Cc1cnc2nccn2c1",      "6-methylimidazo[1,2-a]pyrimidine"),
    ("Cc1ccn2ccnc2n1",      "7-methylimidazo[1,2-a]pyrimidine"),
    # imidazo[1,2-a]pyrazine (CH at 2,3,5,6,8)
    ("c1cn2ccnc2cn1",       "imidazo[1,2-a]pyrazine"),
    ("Cc1cn2ccncc2n1",      "2-methylimidazo[1,2-a]pyrazine"),
    ("Cc1cnc2cnccn12",      "3-methylimidazo[1,2-a]pyrazine"),
    ("Cc1cncc2nccn12",      "5-methylimidazo[1,2-a]pyrazine"),
    ("Cc1cn2ccnc2cn1",      "6-methylimidazo[1,2-a]pyrazine"),
    ("Cc1nccn2ccnc12",      "8-methylimidazo[1,2-a]pyrazine"),
    # imidazo[1,2-b]pyridazine (CH at 2,3,6,7,8)
    ("c1cnn2ccnc2c1",       "imidazo[1,2-b]pyridazine"),
    ("Cc1cn2ncccc2n1",      "2-methylimidazo[1,2-b]pyridazine"),
    ("Cc1cnc2cccnn12",      "3-methylimidazo[1,2-b]pyridazine"),
    ("Cc1ccc2nccn2n1",      "6-methylimidazo[1,2-b]pyridazine"),
    ("Cc1cnn2ccnc2c1",      "7-methylimidazo[1,2-b]pyridazine"),
    ("Cc1ccnn2ccnc12",      "8-methylimidazo[1,2-b]pyridazine"),
    # imidazo[1,2-c]pyrimidine (CH at 2,3,5,7,8)
    ("c1cc2nccn2cn1",       "imidazo[1,2-c]pyrimidine"),
    ("Cc1cn2cnccc2n1",      "2-methylimidazo[1,2-c]pyrimidine"),
    ("Cc1cnc2ccncn12",      "3-methylimidazo[1,2-c]pyrimidine"),
    ("Cc1nccc2nccn12",      "5-methylimidazo[1,2-c]pyrimidine"),
    ("Cc1cc2nccn2cn1",      "7-methylimidazo[1,2-c]pyrimidine"),
    ("Cc1cncn2ccnc12",      "8-methylimidazo[1,2-c]pyrimidine"),
    # imidazo[1,5-a]pyridine (CH at 1,3,5,6,7,8)
    ("c1ccn2cncc2c1",       "imidazo[1,5-a]pyridine"),
    ("Cc1ncn2ccccc12",      "1-methylimidazo[1,5-a]pyridine"),
    ("Cc1ncc2ccccn12",      "3-methylimidazo[1,5-a]pyridine"),
    ("Cc1cccc2cncn12",      "5-methylimidazo[1,5-a]pyridine"),
    ("Cc1ccc2cncn2c1",      "6-methylimidazo[1,5-a]pyridine"),
    ("Cc1ccn2cncc2c1",      "7-methylimidazo[1,5-a]pyridine"),
    ("Cc1cccn2cncc12",      "8-methylimidazo[1,5-a]pyridine"),
    # imidazo[1,5-a]pyrimidine (CH at 2,3,6,8)
    ("c1cnc2cncn2c1",       "imidazo[1,5-a]pyrimidine"),
    ("Cc1ccn2cncc2n1",      "2-methylimidazo[1,5-a]pyrimidine"),
    ("Cc1cnc2cncn2c1",      "3-methylimidazo[1,5-a]pyrimidine"),
    ("Cc1ncc2ncccn12",      "6-methylimidazo[1,5-a]pyrimidine"),
    ("Cc1ncn2cccnc12",      "8-methylimidazo[1,5-a]pyrimidine"),
    # imidazo[1,5-a]pyrazine (CH at 1,3,5,6,8)
    ("c1cn2cncc2cn1",       "imidazo[1,5-a]pyrazine"),
    ("Cc1ncn2ccncc12",      "1-methylimidazo[1,5-a]pyrazine"),
    ("Cc1ncc2cnccn12",      "3-methylimidazo[1,5-a]pyrazine"),
    ("Cc1cncc2cncn12",      "5-methylimidazo[1,5-a]pyrazine"),
    ("Cc1cn2cncc2cn1",      "6-methylimidazo[1,5-a]pyrazine"),
    ("Cc1nccn2cncc12",      "8-methylimidazo[1,5-a]pyrazine"),
    # imidazo[1,5-b]pyridazine (CH at 2,3,5,7)
    ("c1cnn2cncc2c1",       "imidazo[1,5-b]pyridazine"),
    ("Cc1ccc2cncn2n1",      "2-methylimidazo[1,5-b]pyridazine"),
    ("Cc1cnn2cncc2c1",      "3-methylimidazo[1,5-b]pyridazine"),
    ("Cc1ncn2ncccc12",      "5-methylimidazo[1,5-b]pyridazine"),
    ("Cc1ncc2cccnn12",      "7-methylimidazo[1,5-b]pyridazine"),
])
def test_phase700(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
