"""Phase 604: 5-membered heteroaromatics fused with [1,2,3]triazine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[2,3-d][1,2,3]triazine
    ('c1cc2cnnnc2s1',  'thieno[2,3-d][1,2,3]triazine'),
    ('Cc1cc2cnnnc2s1',  '6-methylthieno[2,3-d][1,2,3]triazine'),
    ('Cc1csc2nnncc12',  '5-methylthieno[2,3-d][1,2,3]triazine'),
    ('Cc1nnnc2sccc12',  '4-methylthieno[2,3-d][1,2,3]triazine'),
    #
    # thieno[3,2-d][1,2,3]triazine
    ('c1cc2nnncc2s1',  'thieno[3,2-d][1,2,3]triazine'),
    ('Cc1cc2nnncc2s1',  '6-methylthieno[3,2-d][1,2,3]triazine'),
    ('Cc1csc2cnnnc12',  '7-methylthieno[3,2-d][1,2,3]triazine'),
    ('Cc1nnnc2ccsc12',  '4-methylthieno[3,2-d][1,2,3]triazine'),
    #
    # thieno[3,4-d][1,2,3]triazine
    ('c1nnnc2cscc12',  'thieno[3,4-d][1,2,3]triazine'),
    ('Cc1nnnc2cscc12',  '4-methylthieno[3,4-d][1,2,3]triazine'),
    ('Cc1scc2cnnnc12',  '7-methylthieno[3,4-d][1,2,3]triazine'),
    ('Cc1scc2nnncc12',  '5-methylthieno[3,4-d][1,2,3]triazine'),
    #
    # furo[2,3-d][1,2,3]triazine
    ('c1cc2cnnnc2o1',  'furo[2,3-d][1,2,3]triazine'),
    ('Cc1cc2cnnnc2o1',  '6-methylfuro[2,3-d][1,2,3]triazine'),
    ('Cc1coc2nnncc12',  '5-methylfuro[2,3-d][1,2,3]triazine'),
    ('Cc1nnnc2occc12',  '4-methylfuro[2,3-d][1,2,3]triazine'),
    #
    # furo[3,2-d][1,2,3]triazine
    ('c1cc2nnncc2o1',  'furo[3,2-d][1,2,3]triazine'),
    ('Cc1cc2nnncc2o1',  '6-methylfuro[3,2-d][1,2,3]triazine'),
    ('Cc1coc2cnnnc12',  '7-methylfuro[3,2-d][1,2,3]triazine'),
    ('Cc1nnnc2ccoc12',  '4-methylfuro[3,2-d][1,2,3]triazine'),
    #
    # furo[3,4-d][1,2,3]triazine
    ('c1nnnc2cocc12',  'furo[3,4-d][1,2,3]triazine'),
    ('Cc1nnnc2cocc12',  '4-methylfuro[3,4-d][1,2,3]triazine'),
    ('Cc1occ2cnnnc12',  '7-methylfuro[3,4-d][1,2,3]triazine'),
    ('Cc1occ2nnncc12',  '5-methylfuro[3,4-d][1,2,3]triazine'),
    #
    # isothiazolo[4,3-d][1,2,3]triazine
    ('c1snc2cnnnc12',  'isothiazolo[4,3-d][1,2,3]triazine'),
    ('Cc1snc2cnnnc12',  '7-methylisothiazolo[4,3-d][1,2,3]triazine'),
    ('Cc1nnnc2csnc12',  '4-methylisothiazolo[4,3-d][1,2,3]triazine'),
    #
    # isothiazolo[3,4-d][1,2,3]triazine
    ('c1nnnc2nscc12',  'isothiazolo[3,4-d][1,2,3]triazine'),
    ('Cc1nnnc2nscc12',  '4-methylisothiazolo[3,4-d][1,2,3]triazine'),
    ('Cc1snc2nnncc12',  '5-methylisothiazolo[3,4-d][1,2,3]triazine'),
    #
    # isoxazolo[4,3-d][1,2,3]triazine
    ('c1onc2cnnnc12',  'isoxazolo[4,3-d][1,2,3]triazine'),
    ('Cc1onc2cnnnc12',  '7-methylisoxazolo[4,3-d][1,2,3]triazine'),
    ('Cc1nnnc2conc12',  '4-methylisoxazolo[4,3-d][1,2,3]triazine'),
    #
    # isoxazolo[3,4-d][1,2,3]triazine
    ('c1nnnc2nocc12',  'isoxazolo[3,4-d][1,2,3]triazine'),
    ('Cc1nnnc2nocc12',  '4-methylisoxazolo[3,4-d][1,2,3]triazine'),
    ('Cc1onc2nnncc12',  '5-methylisoxazolo[3,4-d][1,2,3]triazine'),
    #
    # thiazolo[5,4-d][1,2,3]triazine
    ('c1nc2cnnnc2s1',  'thiazolo[5,4-d][1,2,3]triazine'),
    ('Cc1nc2cnnnc2s1',  '6-methylthiazolo[5,4-d][1,2,3]triazine'),
    ('Cc1nnnc2scnc12',  '4-methylthiazolo[5,4-d][1,2,3]triazine'),
    #
    # thiazolo[4,5-d][1,2,3]triazine
    ('c1nc2nnncc2s1',  'thiazolo[4,5-d][1,2,3]triazine'),
    ('Cc1nc2nnncc2s1',  '6-methylthiazolo[4,5-d][1,2,3]triazine'),
    ('Cc1nnnc2ncsc12',  '4-methylthiazolo[4,5-d][1,2,3]triazine'),
    #
    # oxazolo[5,4-d][1,2,3]triazine
    ('c1nc2cnnnc2o1',  'oxazolo[5,4-d][1,2,3]triazine'),
    ('Cc1nc2cnnnc2o1',  '6-methyloxazolo[5,4-d][1,2,3]triazine'),
    ('Cc1nnnc2ocnc12',  '4-methyloxazolo[5,4-d][1,2,3]triazine'),
    #
    # oxazolo[4,5-d][1,2,3]triazine
    ('c1nc2nnncc2o1',  'oxazolo[4,5-d][1,2,3]triazine'),
    ('Cc1nc2nnncc2o1',  '6-methyloxazolo[4,5-d][1,2,3]triazine'),
    ('Cc1nnnc2ncoc12',  '4-methyloxazolo[4,5-d][1,2,3]triazine'),
    #
    # isothiazolo[5,4-d][1,2,3]triazine
    ('c1nnnc2sncc12',  'isothiazolo[5,4-d][1,2,3]triazine'),
    ('Cc1nnnc2sncc12',  '4-methylisothiazolo[5,4-d][1,2,3]triazine'),
    ('Cc1nsc2nnncc12',  '5-methylisothiazolo[5,4-d][1,2,3]triazine'),
    #
    # isoxazolo[4,5-d][1,2,3]triazine
    ('c1noc2cnnnc12',  'isoxazolo[4,5-d][1,2,3]triazine'),
    ('Cc1noc2cnnnc12',  '7-methylisoxazolo[4,5-d][1,2,3]triazine'),
    ('Cc1nnnc2cnoc12',  '4-methylisoxazolo[4,5-d][1,2,3]triazine'),
    #
    # isothiazolo[4,5-d][1,2,3]triazine
    ('c1nsc2cnnnc12',  'isothiazolo[4,5-d][1,2,3]triazine'),
    ('Cc1nsc2cnnnc12',  '7-methylisothiazolo[4,5-d][1,2,3]triazine'),
    ('Cc1nnnc2cnsc12',  '4-methylisothiazolo[4,5-d][1,2,3]triazine'),
    #
    # isoxazolo[5,4-d][1,2,3]triazine
    ('c1nnnc2oncc12',  'isoxazolo[5,4-d][1,2,3]triazine'),
    ('Cc1nnnc2oncc12',  '4-methylisoxazolo[5,4-d][1,2,3]triazine'),
    ('Cc1noc2nnncc12',  '5-methylisoxazolo[5,4-d][1,2,3]triazine'),
])
def test_phase604(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
