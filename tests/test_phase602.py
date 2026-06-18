"""Phase 602: thieno/furo fused with pyridazine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[2,3-c]pyridazine
    ('c1cc2ccsc2nn1',  'thieno[2,3-c]pyridazine'),
    ('Cc1cc2ccsc2nn1',  '3-methylthieno[2,3-c]pyridazine'),
    ('Cc1cnnc2sccc12',  '4-methylthieno[2,3-c]pyridazine'),
    ('Cc1csc2nnccc12',  '5-methylthieno[2,3-c]pyridazine'),
    ('Cc1cc2ccnnc2s1',  '6-methylthieno[2,3-c]pyridazine'),
    #
    # furo[2,3-c]pyridazine
    ('c1cc2ccoc2nn1',  'furo[2,3-c]pyridazine'),
    ('Cc1cc2ccoc2nn1',  '3-methylfuro[2,3-c]pyridazine'),
    ('Cc1cnnc2occc12',  '4-methylfuro[2,3-c]pyridazine'),
    ('Cc1coc2nnccc12',  '5-methylfuro[2,3-c]pyridazine'),
    ('Cc1cc2ccnnc2o1',  '6-methylfuro[2,3-c]pyridazine'),
    #
    # thieno[3,4-c]pyridazine
    ('c1cc2cscc2nn1',  'thieno[3,4-c]pyridazine'),
    ('Cc1cc2cscc2nn1',  '3-methylthieno[3,4-c]pyridazine'),
    ('Cc1cnnc2cscc12',  '4-methylthieno[3,4-c]pyridazine'),
    ('Cc1scc2nnccc12',  '5-methylthieno[3,4-c]pyridazine'),
    ('Cc1scc2ccnnc12',  '7-methylthieno[3,4-c]pyridazine'),
    #
    # furo[3,4-c]pyridazine
    ('c1cc2cocc2nn1',  'furo[3,4-c]pyridazine'),
    ('Cc1cc2cocc2nn1',  '3-methylfuro[3,4-c]pyridazine'),
    ('Cc1cnnc2cocc12',  '4-methylfuro[3,4-c]pyridazine'),
    ('Cc1occ2nnccc12',  '5-methylfuro[3,4-c]pyridazine'),
    ('Cc1occ2ccnnc12',  '7-methylfuro[3,4-c]pyridazine'),
    #
    # thieno[3,2-d]pyridazine
    ('c1cc2cnncc2s1',  'thieno[3,2-d]pyridazine'),
    ('Cc1cc2cnncc2s1',  '2-methylthieno[3,2-d]pyridazine'),
    ('Cc1csc2cnncc12',  '3-methylthieno[3,2-d]pyridazine'),
    ('Cc1nncc2sccc12',  '4-methylthieno[3,2-d]pyridazine'),
    ('Cc1nncc2ccsc12',  '7-methylthieno[3,2-d]pyridazine'),
    #
    # thieno[3,4-d]pyridazine
    ('c1nncc2cscc12',  'thieno[3,4-d]pyridazine'),
    ('Cc1nncc2cscc12',  '1-methylthieno[3,4-d]pyridazine'),
    ('Cc1scc2cnncc12',  '5-methylthieno[3,4-d]pyridazine'),
    #
    # furo[3,2-d]pyridazine
    ('c1cc2cnncc2o1',  'furo[3,2-d]pyridazine'),
    ('Cc1cc2cnncc2o1',  '2-methylfuro[3,2-d]pyridazine'),
    ('Cc1coc2cnncc12',  '3-methylfuro[3,2-d]pyridazine'),
    ('Cc1nncc2occc12',  '4-methylfuro[3,2-d]pyridazine'),
    ('Cc1nncc2ccoc12',  '7-methylfuro[3,2-d]pyridazine'),
    #
    # furo[3,4-d]pyridazine
    ('c1nncc2cocc12',  'furo[3,4-d]pyridazine'),
    ('Cc1nncc2cocc12',  '1-methylfuro[3,4-d]pyridazine'),
    ('Cc1occ2cnncc12',  '5-methylfuro[3,4-d]pyridazine'),
])
def test_phase602(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
