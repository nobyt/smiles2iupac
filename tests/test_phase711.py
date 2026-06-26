"""Phase 711: methyl derivatives of [1,2,5]oxadiazolo and [1,2,5]thiadiazolo series
(8 parent compounds: pyridine, pyrazine, pyridazine, [1,2,4]triazine partners).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # [1,2,5]oxadiazolo[3,4-b]pyridine (CH at 5,6,7)
    ("c1cnc2nonc2c1",   "[1,2,5]oxadiazolo[3,4-b]pyridine"),
    ("Cc1ccc2nonc2n1",  "5-methyl[1,2,5]oxadiazolo[3,4-b]pyridine"),
    ("Cc1cnc2nonc2c1",  "6-methyl[1,2,5]oxadiazolo[3,4-b]pyridine"),
    ("Cc1ccnc2nonc12",  "7-methyl[1,2,5]oxadiazolo[3,4-b]pyridine"),
    # [1,2,5]oxadiazolo[3,4-e]pyrazine (CH at 5)
    ("c1cnc2nonc2n1",   "[1,2,5]oxadiazolo[3,4-e]pyrazine"),
    ("Cc1cnc2nonc2n1",  "5-methyl[1,2,5]oxadiazolo[3,4-e]pyrazine"),
    # [1,2,5]oxadiazolo[3,4-d]pyridazine (CH at 4)
    ("c1nncc2nonc12",   "[1,2,5]oxadiazolo[3,4-d]pyridazine"),
    ("Cc1nncc2nonc12",  "4-methyl[1,2,5]oxadiazolo[3,4-d]pyridazine"),
    # [1,2,5]oxadiazolo[3,4-e][1,2,4]triazine (CH at 6)
    ("c1nnc2nonc2n1",   "[1,2,5]oxadiazolo[3,4-e][1,2,4]triazine"),
    ("Cc1nnc2nonc2n1",  "6-methyl[1,2,5]oxadiazolo[3,4-e][1,2,4]triazine"),
    # [1,2,5]thiadiazolo[3,4-b]pyridine (CH at 5,6,7)
    ("c1cnc2nsnc2c1",   "[1,2,5]thiadiazolo[3,4-b]pyridine"),
    ("Cc1ccc2nsnc2n1",  "5-methyl[1,2,5]thiadiazolo[3,4-b]pyridine"),
    ("Cc1cnc2nsnc2c1",  "6-methyl[1,2,5]thiadiazolo[3,4-b]pyridine"),
    ("Cc1ccnc2nsnc12",  "7-methyl[1,2,5]thiadiazolo[3,4-b]pyridine"),
    # [1,2,5]thiadiazolo[3,4-e]pyrazine (CH at 5)
    ("c1cnc2nsnc2n1",   "[1,2,5]thiadiazolo[3,4-e]pyrazine"),
    ("Cc1cnc2nsnc2n1",  "5-methyl[1,2,5]thiadiazolo[3,4-e]pyrazine"),
    # [1,2,5]thiadiazolo[3,4-d]pyridazine (CH at 4)
    ("c1nncc2nsnc12",   "[1,2,5]thiadiazolo[3,4-d]pyridazine"),
    ("Cc1nncc2nsnc12",  "4-methyl[1,2,5]thiadiazolo[3,4-d]pyridazine"),
    # [1,2,5]thiadiazolo[3,4-e][1,2,4]triazine (CH at 6)
    ("c1nnc2nsnc2n1",   "[1,2,5]thiadiazolo[3,4-e][1,2,4]triazine"),
    ("Cc1nnc2nsnc2n1",  "6-methyl[1,2,5]thiadiazolo[3,4-e][1,2,4]triazine"),
])
def test_phase711(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
