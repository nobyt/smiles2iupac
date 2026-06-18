"""Phase 608: thieno/furo fused with pyridine (b/c-fusion); benzo[c]thiophene; benzo-selenophene."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[3,2-c]pyridine
    ('c1cc2sccc2cn1',  'thieno[3,2-c]pyridine'),
    ('Cc1cc2sccc2cn1',  '6-methylthieno[3,2-c]pyridine'),
    ('Cc1cncc2ccsc12',  '7-methylthieno[3,2-c]pyridine'),
    ('Cc1cc2cnccc2s1',  '2-methylthieno[3,2-c]pyridine'),
    ('Cc1csc2ccncc12',  '3-methylthieno[3,2-c]pyridine'),
    ('Cc1nccc2sccc12',  '4-methylthieno[3,2-c]pyridine'),
    #
    # furo[3,2-b]pyridine
    ('c1cnc2ccoc2c1',  'furo[3,2-b]pyridine'),
    ('Cc1cnc2ccoc2c1',  '6-methylfuro[3,2-b]pyridine'),
    ('Cc1ccc2occc2n1',  '5-methylfuro[3,2-b]pyridine'),
    ('Cc1coc2cccnc12',  '3-methylfuro[3,2-b]pyridine'),
    ('Cc1cc2ncccc2o1',  '2-methylfuro[3,2-b]pyridine'),
    ('Cc1ccnc2ccoc12',  '7-methylfuro[3,2-b]pyridine'),
    #
    # furo[3,4-b]pyridine
    ('c1cnc2cocc2c1',  'furo[3,4-b]pyridine'),
    ('Cc1cnc2cocc2c1',  '3-methylfuro[3,4-b]pyridine'),
    ('Cc1ccc2cocc2n1',  '2-methylfuro[3,4-b]pyridine'),
    ('Cc1occ2cccnc12',  '7-methylfuro[3,4-b]pyridine'),
    ('Cc1occ2ncccc12',  '5-methylfuro[3,4-b]pyridine'),
    ('Cc1ccnc2cocc12',  '4-methylfuro[3,4-b]pyridine'),
    #
    # benzo[c]thiophene
    ('c1ccc2cscc2c1',  'benzo[c]thiophene'),
    ('Cc1ccc2cscc2c1',  '5-methylbenzo[c]thiophene'),
    ('Cc1cccc2cscc12',  '4-methylbenzo[c]thiophene'),
    ('Cc1scc2ccccc12',  '1-methylbenzo[c]thiophene'),
    #
    # benzo[b]selenophene
    ('c1ccc2[se]ccc2c1',  'benzo[b]selenophene'),
    ('Cc1ccc2[se]ccc2c1',  '5-methylbenzo[b]selenophene'),
    ('Cc1ccc2cc[se]c2c1',  '6-methylbenzo[b]selenophene'),
    ('Cc1cccc2cc[se]c12',  '7-methylbenzo[b]selenophene'),
    ('Cc1cc2ccccc2[se]1',  '2-methylbenzo[b]selenophene'),
    ('Cc1c[se]c2ccccc12',  '3-methylbenzo[b]selenophene'),
    ('Cc1cccc2[se]ccc12',  '4-methylbenzo[b]selenophene'),
    #
    # benzo[c]selenophene
    ('c1ccc2c[se]cc2c1',  'benzo[c]selenophene'),
    ('Cc1ccc2c[se]cc2c1',  '5-methylbenzo[c]selenophene'),
    ('Cc1cccc2c[se]cc12',  '4-methylbenzo[c]selenophene'),
    ('Cc1[se]cc2ccccc12',  '1-methylbenzo[c]selenophene'),
    #
    # furo[3,2-c]pyridine
    ('c1cc2occc2cn1',  'furo[3,2-c]pyridine'),
    ('Cc1cc2occc2cn1',  '6-methylfuro[3,2-c]pyridine'),
    ('Cc1cncc2ccoc12',  '7-methylfuro[3,2-c]pyridine'),
    ('Cc1cc2cnccc2o1',  '2-methylfuro[3,2-c]pyridine'),
    ('Cc1coc2ccncc12',  '3-methylfuro[3,2-c]pyridine'),
    ('Cc1nccc2occc12',  '4-methylfuro[3,2-c]pyridine'),
    #
    # furo[2,3-e]pyridazine
    ('c1cc2occc2nn1',  'furo[2,3-e]pyridazine'),
    ('Cc1cc2occc2nn1',  '3-methylfuro[2,3-e]pyridazine'),
    ('Cc1cnnc2ccoc12',  '4-methylfuro[2,3-e]pyridazine'),
    ('Cc1cc2nnccc2o1',  '6-methylfuro[2,3-e]pyridazine'),
    ('Cc1coc2ccnnc12',  '7-methylfuro[2,3-e]pyridazine'),
    #
    # thieno[2,3-c]pyridine
    ('c1cc2ccsc2cn1',  'thieno[2,3-c]pyridine'),
    ('Cc1cc2ccsc2cn1',  '5-methylthieno[2,3-c]pyridine'),
    ('Cc1cncc2sccc12',  '4-methylthieno[2,3-c]pyridine'),
    ('Cc1csc2cnccc12',  '3-methylthieno[2,3-c]pyridine'),
    ('Cc1cc2ccncc2s1',  '2-methylthieno[2,3-c]pyridine'),
    ('Cc1nccc2ccsc12',  '7-methylthieno[2,3-c]pyridine'),
    #
    # furo[2,3-c]pyridine
    ('c1cc2ccoc2cn1',  'furo[2,3-c]pyridine'),
    ('Cc1cc2ccoc2cn1',  '5-methylfuro[2,3-c]pyridine'),
    ('Cc1cncc2occc12',  '4-methylfuro[2,3-c]pyridine'),
    ('Cc1coc2cnccc12',  '3-methylfuro[2,3-c]pyridine'),
    ('Cc1cc2ccncc2o1',  '2-methylfuro[2,3-c]pyridine'),
    ('Cc1nccc2ccoc12',  '7-methylfuro[2,3-c]pyridine'),
    #
    # furo[3,4-c]pyridine
    ('c1cc2cocc2cn1',  'furo[3,4-c]pyridine'),
    ('Cc1cc2cocc2cn1',  '6-methylfuro[3,4-c]pyridine'),
    ('Cc1cncc2cocc12',  '7-methylfuro[3,4-c]pyridine'),
    ('Cc1occ2cnccc12',  '1-methylfuro[3,4-c]pyridine'),
    ('Cc1occ2ccncc12',  '3-methylfuro[3,4-c]pyridine'),
    ('Cc1nccc2cocc12',  '4-methylfuro[3,4-c]pyridine'),
    #
    # thieno[3,2-c]pyridazine
    ('c1cc2sccc2nn1',  'thieno[3,2-c]pyridazine'),
    ('Cc1cc2sccc2nn1',  '3-methylthieno[3,2-c]pyridazine'),
    ('Cc1cnnc2ccsc12',  '4-methylthieno[3,2-c]pyridazine'),
    ('Cc1cc2nnccc2s1',  '6-methylthieno[3,2-c]pyridazine'),
    ('Cc1csc2ccnnc12',  '7-methylthieno[3,2-c]pyridazine'),
])
def test_phase608(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
