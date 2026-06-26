"""Phase 715: methyl derivatives of oxazolo/thiazolo[x,y-b/c] and [x,y-d][1,2,3]triazine series
(16 parent compounds: 4,5- and 5,4- b/c-fused with pyridine/pyridazine
and [1,2,3]triazine partners for both oxazolo and thiazolo).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # oxazolo[4,5-b]pyridine (CH at 2,5,6,7)
    ("c1cnc2ncoc2c1",       "oxazolo[4,5-b]pyridine"),
    ("Cc1nc2ncccc2o1",      "2-methyloxazolo[4,5-b]pyridine"),
    ("Cc1ccc2ocnc2n1",      "5-methyloxazolo[4,5-b]pyridine"),
    ("Cc1cnc2ncoc2c1",      "6-methyloxazolo[4,5-b]pyridine"),
    ("Cc1ccnc2ncoc12",      "7-methyloxazolo[4,5-b]pyridine"),
    # oxazolo[4,5-c]pyridine (CH at 2,4,6,7)
    ("c1cc2ocnc2cn1",       "oxazolo[4,5-c]pyridine"),
    ("Cc1nc2cnccc2o1",      "2-methyloxazolo[4,5-c]pyridine"),
    ("Cc1nccc2ocnc12",      "4-methyloxazolo[4,5-c]pyridine"),
    ("Cc1cc2ocnc2cn1",      "6-methyloxazolo[4,5-c]pyridine"),
    ("Cc1cncc2ncoc12",      "7-methyloxazolo[4,5-c]pyridine"),
    # oxazolo[4,5-c]pyridazine (CH at 3,4,6)
    ("c1cc2ocnc2nn1",       "oxazolo[4,5-c]pyridazine"),
    ("Cc1cc2ocnc2nn1",      "3-methyloxazolo[4,5-c]pyridazine"),
    ("Cc1cnnc2ncoc12",      "4-methyloxazolo[4,5-c]pyridazine"),
    ("Cc1nc2nnccc2o1",      "6-methyloxazolo[4,5-c]pyridazine"),
    # oxazolo[4,5-d][1,2,3]triazine (CH at 4,6)
    ("c1nc2nnncc2o1",       "oxazolo[4,5-d][1,2,3]triazine"),
    ("Cc1nnnc2ncoc12",      "4-methyloxazolo[4,5-d][1,2,3]triazine"),
    ("Cc1nc2nnncc2o1",      "6-methyloxazolo[4,5-d][1,2,3]triazine"),
    # oxazolo[5,4-b]pyridine (CH at 2,5,6,7)
    ("c1cnc2ocnc2c1",       "oxazolo[5,4-b]pyridine"),
    ("Cc1nc2cccnc2o1",      "2-methyloxazolo[5,4-b]pyridine"),
    ("Cc1ccc2ncoc2n1",      "5-methyloxazolo[5,4-b]pyridine"),
    ("Cc1cnc2ocnc2c1",      "6-methyloxazolo[5,4-b]pyridine"),
    ("Cc1ccnc2ocnc12",      "7-methyloxazolo[5,4-b]pyridine"),
    # oxazolo[5,4-c]pyridine (CH at 2,4,6,7)
    ("c1cc2ncoc2cn1",       "oxazolo[5,4-c]pyridine"),
    ("Cc1nc2ccncc2o1",      "2-methyloxazolo[5,4-c]pyridine"),
    ("Cc1nccc2ncoc12",      "4-methyloxazolo[5,4-c]pyridine"),
    ("Cc1cc2ncoc2cn1",      "6-methyloxazolo[5,4-c]pyridine"),
    ("Cc1cncc2ocnc12",      "7-methyloxazolo[5,4-c]pyridine"),
    # oxazolo[5,4-c]pyridazine (CH at 3,4,6)
    ("c1cc2ncoc2nn1",       "oxazolo[5,4-c]pyridazine"),
    ("Cc1cc2ncoc2nn1",      "3-methyloxazolo[5,4-c]pyridazine"),
    ("Cc1cnnc2ocnc12",      "4-methyloxazolo[5,4-c]pyridazine"),
    ("Cc1nc2ccnnc2o1",      "6-methyloxazolo[5,4-c]pyridazine"),
    # oxazolo[5,4-d][1,2,3]triazine (CH at 4,6)
    ("c1nc2cnnnc2o1",       "oxazolo[5,4-d][1,2,3]triazine"),
    ("Cc1nnnc2ocnc12",      "4-methyloxazolo[5,4-d][1,2,3]triazine"),
    ("Cc1nc2cnnnc2o1",      "6-methyloxazolo[5,4-d][1,2,3]triazine"),
    # thiazolo[4,5-b]pyridine (CH at 2,5,6,7)
    ("c1cnc2ncsc2c1",       "thiazolo[4,5-b]pyridine"),
    ("Cc1nc2ncccc2s1",      "2-methylthiazolo[4,5-b]pyridine"),
    ("Cc1ccc2scnc2n1",      "5-methylthiazolo[4,5-b]pyridine"),
    ("Cc1cnc2ncsc2c1",      "6-methylthiazolo[4,5-b]pyridine"),
    ("Cc1ccnc2ncsc12",      "7-methylthiazolo[4,5-b]pyridine"),
    # thiazolo[4,5-c]pyridine (CH at 2,4,6,7)
    ("c1cc2scnc2cn1",       "thiazolo[4,5-c]pyridine"),
    ("Cc1nc2cnccc2s1",      "2-methylthiazolo[4,5-c]pyridine"),
    ("Cc1nccc2scnc12",      "4-methylthiazolo[4,5-c]pyridine"),
    ("Cc1cc2scnc2cn1",      "6-methylthiazolo[4,5-c]pyridine"),
    ("Cc1cncc2ncsc12",      "7-methylthiazolo[4,5-c]pyridine"),
    # thiazolo[4,5-c]pyridazine (CH at 3,4,6)
    ("c1cc2scnc2nn1",       "thiazolo[4,5-c]pyridazine"),
    ("Cc1cc2scnc2nn1",      "3-methylthiazolo[4,5-c]pyridazine"),
    ("Cc1cnnc2ncsc12",      "4-methylthiazolo[4,5-c]pyridazine"),
    ("Cc1nc2nnccc2s1",      "6-methylthiazolo[4,5-c]pyridazine"),
    # thiazolo[4,5-d][1,2,3]triazine (CH at 4,6)
    ("c1nc2nnncc2s1",       "thiazolo[4,5-d][1,2,3]triazine"),
    ("Cc1nnnc2ncsc12",      "4-methylthiazolo[4,5-d][1,2,3]triazine"),
    ("Cc1nc2nnncc2s1",      "6-methylthiazolo[4,5-d][1,2,3]triazine"),
    # thiazolo[5,4-b]pyridine (CH at 2,5,6,7)
    ("c1cnc2scnc2c1",       "thiazolo[5,4-b]pyridine"),
    ("Cc1nc2cccnc2s1",      "2-methylthiazolo[5,4-b]pyridine"),
    ("Cc1ccc2ncsc2n1",      "5-methylthiazolo[5,4-b]pyridine"),
    ("Cc1cnc2scnc2c1",      "6-methylthiazolo[5,4-b]pyridine"),
    ("Cc1ccnc2scnc12",      "7-methylthiazolo[5,4-b]pyridine"),
    # thiazolo[5,4-c]pyridine (CH at 2,4,6,7)
    ("c1cc2ncsc2cn1",       "thiazolo[5,4-c]pyridine"),
    ("Cc1nc2ccncc2s1",      "2-methylthiazolo[5,4-c]pyridine"),
    ("Cc1nccc2ncsc12",      "4-methylthiazolo[5,4-c]pyridine"),
    ("Cc1cc2ncsc2cn1",      "6-methylthiazolo[5,4-c]pyridine"),
    ("Cc1cncc2scnc12",      "7-methylthiazolo[5,4-c]pyridine"),
    # thiazolo[5,4-c]pyridazine (CH at 3,4,6)
    ("c1cc2ncsc2nn1",       "thiazolo[5,4-c]pyridazine"),
    ("Cc1cc2ncsc2nn1",      "3-methylthiazolo[5,4-c]pyridazine"),
    ("Cc1cnnc2scnc12",      "4-methylthiazolo[5,4-c]pyridazine"),
    ("Cc1nc2ccnnc2s1",      "6-methylthiazolo[5,4-c]pyridazine"),
    # thiazolo[5,4-d][1,2,3]triazine (CH at 4,6)
    ("c1nc2cnnnc2s1",       "thiazolo[5,4-d][1,2,3]triazine"),
    ("Cc1nnnc2scnc12",      "4-methylthiazolo[5,4-d][1,2,3]triazine"),
    ("Cc1nc2cnnnc2s1",      "6-methylthiazolo[5,4-d][1,2,3]triazine"),
])
def test_phase715(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
