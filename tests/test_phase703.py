"""Phase 703: methyl derivatives of [1,2,3]triazolo[1,5-a/b] and tetrazolo[1,5-a/b/d] series
(9 parent compounds).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # [1,2,3]triazolo[1,5-a]pyridine (CH at 3,5,6,7)
    ("c1ccn2nncc2c1",      "[1,2,3]triazolo[1,5-a]pyridine"),
    ("Cc1nnn2ccccc12",     "3-methyl[1,2,3]triazolo[1,5-a]pyridine"),
    ("Cc1ccn2nncc2c1",     "5-methyl[1,2,3]triazolo[1,5-a]pyridine"),
    ("Cc1ccc2cnnn2c1",     "6-methyl[1,2,3]triazolo[1,5-a]pyridine"),
    ("Cc1cccc2cnnn12",     "7-methyl[1,2,3]triazolo[1,5-a]pyridine"),
    # [1,2,3]triazolo[1,5-a]pyrimidine (CH at 3,5,6,7)
    ("c1cnc2cnnn2c1",      "[1,2,3]triazolo[1,5-a]pyrimidine"),
    ("Cc1nnn2cccnc12",     "3-methyl[1,2,3]triazolo[1,5-a]pyrimidine"),
    ("Cc1ccn2nncc2n1",     "5-methyl[1,2,3]triazolo[1,5-a]pyrimidine"),
    ("Cc1cnc2cnnn2c1",     "6-methyl[1,2,3]triazolo[1,5-a]pyrimidine"),
    ("Cc1ccnc2cnnn12",     "7-methyl[1,2,3]triazolo[1,5-a]pyrimidine"),
    # [1,2,3]triazolo[1,5-b]pyridazine (CH at 3,5,6)
    ("c1cnn2nncc2c1",      "[1,2,3]triazolo[1,5-b]pyridazine"),
    ("Cc1nnn2ncccc12",     "3-methyl[1,2,3]triazolo[1,5-b]pyridazine"),
    ("Cc1cnn2nncc2c1",     "5-methyl[1,2,3]triazolo[1,5-b]pyridazine"),
    ("Cc1ccc2cnnn2n1",     "6-methyl[1,2,3]triazolo[1,5-b]pyridazine"),
    # tetrazolo[1,5-a]pyridine (CH at 5,6,7,8)
    ("c1ccn2nnnc2c1",      "tetrazolo[1,5-a]pyridine"),
    ("Cc1cccc2nnnn12",     "5-methyltetrazolo[1,5-a]pyridine"),
    ("Cc1ccc2nnnn2c1",     "6-methyltetrazolo[1,5-a]pyridine"),
    ("Cc1ccn2nnnc2c1",     "7-methyltetrazolo[1,5-a]pyridine"),
    ("Cc1cccn2nnnc12",     "8-methyltetrazolo[1,5-a]pyridine"),
    # tetrazolo[1,5-a]pyrimidine (CH at 5,6,7)
    ("c1cnc2nnnn2c1",      "tetrazolo[1,5-a]pyrimidine"),
    ("Cc1ccn2nnnc2n1",     "5-methyltetrazolo[1,5-a]pyrimidine"),
    ("Cc1cnc2nnnn2c1",     "6-methyltetrazolo[1,5-a]pyrimidine"),
    ("Cc1ccnc2nnnn12",     "7-methyltetrazolo[1,5-a]pyrimidine"),
    # tetrazolo[1,5-a]pyrazine (CH at 5,6,8)
    ("c1cn2nnnc2cn1",      "tetrazolo[1,5-a]pyrazine"),
    ("Cc1cncc2nnnn12",     "5-methyltetrazolo[1,5-a]pyrazine"),
    ("Cc1cn2nnnc2cn1",     "6-methyltetrazolo[1,5-a]pyrazine"),
    ("Cc1nccn2nnnc12",     "8-methyltetrazolo[1,5-a]pyrazine"),
    # tetrazolo[1,5-b]pyridazine (CH at 6,7,8)
    ("c1cnn2nnnc2c1",      "tetrazolo[1,5-b]pyridazine"),
    ("Cc1ccc2nnnn2n1",     "6-methyltetrazolo[1,5-b]pyridazine"),
    ("Cc1cnn2nnnc2c1",     "7-methyltetrazolo[1,5-b]pyridazine"),
    ("Cc1ccnn2nnnc12",     "8-methyltetrazolo[1,5-b]pyridazine"),
    # tetrazolo[1,5-b][1,2,4]triazine (CH at 6,7)
    ("c1cnn2nnnc2n1",      "tetrazolo[1,5-b][1,2,4]triazine"),
    ("Cc1cnc2nnnn2n1",     "6-methyltetrazolo[1,5-b][1,2,4]triazine"),
    ("Cc1cnn2nnnc2n1",     "7-methyltetrazolo[1,5-b][1,2,4]triazine"),
    # tetrazolo[1,5-d][1,2,4]triazine (CH at 5,8)
    ("c1nncn2nnnc12",      "tetrazolo[1,5-d][1,2,4]triazine"),
    ("Cc1nncc2nnnn12",     "5-methyltetrazolo[1,5-d][1,2,4]triazine"),
    ("Cc1nncn2nnnc12",     "8-methyltetrazolo[1,5-d][1,2,4]triazine"),
])
def test_phase703(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
