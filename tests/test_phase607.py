"""Phase 607: tetrazolo fusions; isothiazolo/isoxazolo/oxazolo/thiazolo-pyridazine; pyrrolo/triazolo-triazine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # tetrazolo[1,5-a]pyrazine
    ('c1cn2nnnc2cn1',  'tetrazolo[1,5-a]pyrazine'),
    ('Cc1cn2nnnc2cn1',  '6-methyltetrazolo[1,5-a]pyrazine'),
    ('Cc1cncc2nnnn12',  '5-methyltetrazolo[1,5-a]pyrazine'),
    ('Cc1nccn2nnnc12',  '8-methyltetrazolo[1,5-a]pyrazine'),
    #
    # tetrazolo[1,5-b][1,2,4]triazine
    ('c1cnn2nnnc2n1',  'tetrazolo[1,5-b][1,2,4]triazine'),
    ('Cc1cnn2nnnc2n1',  '7-methyltetrazolo[1,5-b][1,2,4]triazine'),
    ('Cc1cnc2nnnn2n1',  '6-methyltetrazolo[1,5-b][1,2,4]triazine'),
    #
    # tetrazolo[1,5-d][1,2,4]triazine
    ('c1nncn2nnnc12',  'tetrazolo[1,5-d][1,2,4]triazine'),
    ('Cc1nncn2nnnc12',  '8-methyltetrazolo[1,5-d][1,2,4]triazine'),
    ('Cc1nncc2nnnn12',  '5-methyltetrazolo[1,5-d][1,2,4]triazine'),
    #
    # isothiazolo[5,4-d]pyridazine
    ('c1nncc2sncc12',  'isothiazolo[5,4-d]pyridazine'),
    ('Cc1nncc2sncc12',  '4-methylisothiazolo[5,4-d]pyridazine'),
    ('Cc1nncc2cnsc12',  '7-methylisothiazolo[5,4-d]pyridazine'),
    ('Cc1nsc2cnncc12',  '3-methylisothiazolo[5,4-d]pyridazine'),
    #
    # isothiazolo[3,4-d]pyridazine
    ('c1nncc2nscc12',  'isothiazolo[3,4-d]pyridazine'),
    ('Cc1nncc2nscc12',  '4-methylisothiazolo[3,4-d]pyridazine'),
    ('Cc1nncc2csnc12',  '7-methylisothiazolo[3,4-d]pyridazine'),
    ('Cc1snc2cnncc12',  '3-methylisothiazolo[3,4-d]pyridazine'),
    #
    # isoxazolo[5,4-d]pyridazine
    ('c1nncc2oncc12',  'isoxazolo[5,4-d]pyridazine'),
    ('Cc1nncc2oncc12',  '4-methylisoxazolo[5,4-d]pyridazine'),
    ('Cc1nncc2cnoc12',  '7-methylisoxazolo[5,4-d]pyridazine'),
    ('Cc1noc2cnncc12',  '3-methylisoxazolo[5,4-d]pyridazine'),
    #
    # isoxazolo[3,4-d]pyridazine
    ('c1nncc2nocc12',  'isoxazolo[3,4-d]pyridazine'),
    ('Cc1nncc2nocc12',  '4-methylisoxazolo[3,4-d]pyridazine'),
    ('Cc1nncc2conc12',  '7-methylisoxazolo[3,4-d]pyridazine'),
    ('Cc1onc2cnncc12',  '3-methylisoxazolo[3,4-d]pyridazine'),
    #
    # oxazolo[4,5-d]pyridazine
    ('c1nc2cnncc2o1',  'oxazolo[4,5-d]pyridazine'),
    ('Cc1nc2cnncc2o1',  '2-methyloxazolo[4,5-d]pyridazine'),
    ('Cc1nncc2ocnc12',  '4-methyloxazolo[4,5-d]pyridazine'),
    ('Cc1nncc2ncoc12',  '7-methyloxazolo[4,5-d]pyridazine'),
    #
    # thiazolo[4,5-d]pyridazine
    ('c1nc2cnncc2s1',  'thiazolo[4,5-d]pyridazine'),
    ('Cc1nc2cnncc2s1',  '2-methylthiazolo[4,5-d]pyridazine'),
    ('Cc1nncc2scnc12',  '4-methylthiazolo[4,5-d]pyridazine'),
    ('Cc1nncc2ncsc12',  '7-methylthiazolo[4,5-d]pyridazine'),
    #
    # [1,2,4]triazolo[1,5-b][1,2,4]triazine
    ('c1cnn2ncnc2n1',  '[1,2,4]triazolo[1,5-b][1,2,4]triazine'),
    ('Cc1cnn2ncnc2n1',  '7-methyl[1,2,4]triazolo[1,5-b][1,2,4]triazine'),
    ('Cc1cnc2ncnn2n1',  '6-methyl[1,2,4]triazolo[1,5-b][1,2,4]triazine'),
    ('Cc1nc2nccnn2n1',  '2-methyl[1,2,4]triazolo[1,5-b][1,2,4]triazine'),
    #
    # pyrrolo[1,2-b][1,2,4]triazine
    ('c1cc2nccnn2c1',  'pyrrolo[1,2-b][1,2,4]triazine'),
    ('Cc1cc2nccnn2c1',  '7-methylpyrrolo[1,2-b][1,2,4]triazine'),
    ('Cc1ccn2nccnc12',  '8-methylpyrrolo[1,2-b][1,2,4]triazine'),
    ('Cc1cnn2cccc2n1',  '2-methylpyrrolo[1,2-b][1,2,4]triazine'),
    ('Cc1cnc2cccn2n1',  '3-methylpyrrolo[1,2-b][1,2,4]triazine'),
    ('Cc1ccc2nccnn12',  '6-methylpyrrolo[1,2-b][1,2,4]triazine'),
    #
    # tetrazolo[1,5-a]pyrimidine
    ('c1cnc2nnnn2c1',  'tetrazolo[1,5-a]pyrimidine'),
    ('Cc1cnc2nnnn2c1',  '6-methyltetrazolo[1,5-a]pyrimidine'),
    ('Cc1ccn2nnnc2n1',  '5-methyltetrazolo[1,5-a]pyrimidine'),
    ('Cc1ccnc2nnnn12',  '7-methyltetrazolo[1,5-a]pyrimidine'),
])
def test_phase607(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
