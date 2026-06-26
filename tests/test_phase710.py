"""Phase 710: methyl derivatives of [1,2,3]oxadiazolo and [1,2,3]thiadiazolo series
(14 parent compounds: pyridine, pyridine-c, pyrazine, pyridazine,
[1,2,4]triazine, [1,2,3]triazine partners).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # [1,2,3]oxadiazolo[4,5-b]pyridine (CH at 5,6,7)
    ("c1cnc2nnoc2c1",   "[1,2,3]oxadiazolo[4,5-b]pyridine"),
    ("Cc1ccc2onnc2n1",  "5-methyl[1,2,3]oxadiazolo[4,5-b]pyridine"),
    ("Cc1cnc2nnoc2c1",  "6-methyl[1,2,3]oxadiazolo[4,5-b]pyridine"),
    ("Cc1ccnc2nnoc12",  "7-methyl[1,2,3]oxadiazolo[4,5-b]pyridine"),
    # [1,2,3]oxadiazolo[4,5-c]pyridine (CH at 4,6,7)
    ("c1cc2onnc2cn1",   "[1,2,3]oxadiazolo[4,5-c]pyridine"),
    ("Cc1nccc2onnc12",  "4-methyl[1,2,3]oxadiazolo[4,5-c]pyridine"),
    ("Cc1cc2onnc2cn1",  "6-methyl[1,2,3]oxadiazolo[4,5-c]pyridine"),
    ("Cc1cncc2nnoc12",  "7-methyl[1,2,3]oxadiazolo[4,5-c]pyridine"),
    # [1,2,3]oxadiazolo[4,5-e]pyrazine (CH at 5,6)
    ("c1cnc2onnc2n1",   "[1,2,3]oxadiazolo[4,5-e]pyrazine"),
    ("Cc1cnc2onnc2n1",  "5-methyl[1,2,3]oxadiazolo[4,5-e]pyrazine"),
    ("Cc1cnc2nnoc2n1",  "6-methyl[1,2,3]oxadiazolo[4,5-e]pyrazine"),
    # [1,2,3]oxadiazolo[5,4-d]pyridazine (CH at 4,7)
    ("c1nncc2onnc12",   "[1,2,3]oxadiazolo[5,4-d]pyridazine"),
    ("Cc1nncc2onnc12",  "4-methyl[1,2,3]oxadiazolo[5,4-d]pyridazine"),
    ("Cc1nncc2nnoc12",  "7-methyl[1,2,3]oxadiazolo[5,4-d]pyridazine"),
    # [1,2,3]oxadiazolo[4,5-e][1,2,4]triazine (CH at 6)
    ("c1nnc2onnc2n1",   "[1,2,3]oxadiazolo[4,5-e][1,2,4]triazine"),
    ("Cc1nnc2onnc2n1",  "6-methyl[1,2,3]oxadiazolo[4,5-e][1,2,4]triazine"),
    # [1,2,3]oxadiazolo[5,4-d][1,2,3]triazine (CH at 7)
    ("c1nnnc2onnc12",   "[1,2,3]oxadiazolo[5,4-d][1,2,3]triazine"),
    ("Cc1nnnc2onnc12",  "7-methyl[1,2,3]oxadiazolo[5,4-d][1,2,3]triazine"),
    # [1,2,3]oxadiazolo[4,5-d][1,2,3]triazine (CH at 7)
    ("c1nnnc2nnoc12",   "[1,2,3]oxadiazolo[4,5-d][1,2,3]triazine"),
    ("Cc1nnnc2nnoc12",  "7-methyl[1,2,3]oxadiazolo[4,5-d][1,2,3]triazine"),
    # [1,2,3]thiadiazolo[4,5-b]pyridine (CH at 5,6,7)
    ("c1cnc2nnsc2c1",   "[1,2,3]thiadiazolo[4,5-b]pyridine"),
    ("Cc1ccc2snnc2n1",  "5-methyl[1,2,3]thiadiazolo[4,5-b]pyridine"),
    ("Cc1cnc2nnsc2c1",  "6-methyl[1,2,3]thiadiazolo[4,5-b]pyridine"),
    ("Cc1ccnc2nnsc12",  "7-methyl[1,2,3]thiadiazolo[4,5-b]pyridine"),
    # [1,2,3]thiadiazolo[4,5-c]pyridine (CH at 4,6,7)
    ("c1cc2snnc2cn1",   "[1,2,3]thiadiazolo[4,5-c]pyridine"),
    ("Cc1nccc2snnc12",  "4-methyl[1,2,3]thiadiazolo[4,5-c]pyridine"),
    ("Cc1cc2snnc2cn1",  "6-methyl[1,2,3]thiadiazolo[4,5-c]pyridine"),
    ("Cc1cncc2nnsc12",  "7-methyl[1,2,3]thiadiazolo[4,5-c]pyridine"),
    # [1,2,3]thiadiazolo[4,5-e]pyrazine (CH at 5,6)
    ("c1cnc2snnc2n1",   "[1,2,3]thiadiazolo[4,5-e]pyrazine"),
    ("Cc1cnc2snnc2n1",  "5-methyl[1,2,3]thiadiazolo[4,5-e]pyrazine"),
    ("Cc1cnc2nnsc2n1",  "6-methyl[1,2,3]thiadiazolo[4,5-e]pyrazine"),
    # [1,2,3]thiadiazolo[5,4-d]pyridazine (CH at 4,7)
    ("c1nncc2snnc12",   "[1,2,3]thiadiazolo[5,4-d]pyridazine"),
    ("Cc1nncc2snnc12",  "4-methyl[1,2,3]thiadiazolo[5,4-d]pyridazine"),
    ("Cc1nncc2nnsc12",  "7-methyl[1,2,3]thiadiazolo[5,4-d]pyridazine"),
    # [1,2,3]thiadiazolo[4,5-e][1,2,4]triazine (CH at 6)
    ("c1nnc2snnc2n1",   "[1,2,3]thiadiazolo[4,5-e][1,2,4]triazine"),
    ("Cc1nnc2snnc2n1",  "6-methyl[1,2,3]thiadiazolo[4,5-e][1,2,4]triazine"),
    # [1,2,3]thiadiazolo[5,4-d][1,2,3]triazine (CH at 7)
    ("c1nnnc2snnc12",   "[1,2,3]thiadiazolo[5,4-d][1,2,3]triazine"),
    ("Cc1nnnc2snnc12",  "7-methyl[1,2,3]thiadiazolo[5,4-d][1,2,3]triazine"),
    # [1,2,3]thiadiazolo[4,5-d][1,2,3]triazine (CH at 7)
    ("c1nnnc2nnsc12",   "[1,2,3]thiadiazolo[4,5-d][1,2,3]triazine"),
    ("Cc1nnnc2nnsc12",  "7-methyl[1,2,3]thiadiazolo[4,5-d][1,2,3]triazine"),
])
def test_phase710(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
