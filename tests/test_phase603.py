"""Phase 603: 5-membered heteroaromatics fused with [1,2,4]triazine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[2,3-e][1,2,4]triazine
    ('c1nnc2ccsc2n1',  'thieno[2,3-e][1,2,4]triazine'),
    ('Cc1nnc2ccsc2n1',  '3-methylthieno[2,3-e][1,2,4]triazine'),
    ('Cc1csc2ncnnc12',  '7-methylthieno[2,3-e][1,2,4]triazine'),
    ('Cc1cc2nncnc2s1',  '6-methylthieno[2,3-e][1,2,4]triazine'),
    #
    # thieno[3,2-e][1,2,4]triazine
    ('c1nnc2sccc2n1',  'thieno[3,2-e][1,2,4]triazine'),
    ('Cc1nnc2sccc2n1',  '3-methylthieno[3,2-e][1,2,4]triazine'),
    ('Cc1cc2ncnnc2s1',  '6-methylthieno[3,2-e][1,2,4]triazine'),
    ('Cc1csc2nncnc12',  '5-methylthieno[3,2-e][1,2,4]triazine'),
    #
    # thieno[3,4-e][1,2,4]triazine
    ('c1nnc2cscc2n1',  'thieno[3,4-e][1,2,4]triazine'),
    ('Cc1nnc2cscc2n1',  '3-methylthieno[3,4-e][1,2,4]triazine'),
    ('Cc1scc2ncnnc12',  '7-methylthieno[3,4-e][1,2,4]triazine'),
    ('Cc1scc2nncnc12',  '5-methylthieno[3,4-e][1,2,4]triazine'),
    #
    # furo[2,3-e][1,2,4]triazine
    ('c1nnc2ccoc2n1',  'furo[2,3-e][1,2,4]triazine'),
    ('Cc1nnc2ccoc2n1',  '3-methylfuro[2,3-e][1,2,4]triazine'),
    ('Cc1coc2ncnnc12',  '7-methylfuro[2,3-e][1,2,4]triazine'),
    ('Cc1cc2nncnc2o1',  '6-methylfuro[2,3-e][1,2,4]triazine'),
    #
    # furo[3,2-e][1,2,4]triazine
    ('c1nnc2occc2n1',  'furo[3,2-e][1,2,4]triazine'),
    ('Cc1nnc2occc2n1',  '3-methylfuro[3,2-e][1,2,4]triazine'),
    ('Cc1cc2ncnnc2o1',  '6-methylfuro[3,2-e][1,2,4]triazine'),
    ('Cc1coc2nncnc12',  '5-methylfuro[3,2-e][1,2,4]triazine'),
    #
    # furo[3,4-e][1,2,4]triazine
    ('c1nnc2cocc2n1',  'furo[3,4-e][1,2,4]triazine'),
    ('Cc1nnc2cocc2n1',  '3-methylfuro[3,4-e][1,2,4]triazine'),
    ('Cc1occ2ncnnc12',  '7-methylfuro[3,4-e][1,2,4]triazine'),
    ('Cc1occ2nncnc12',  '5-methylfuro[3,4-e][1,2,4]triazine'),
    #
    # isothiazolo[3,4-e][1,2,4]triazine
    ('c1nnc2csnc2n1',  'isothiazolo[3,4-e][1,2,4]triazine'),
    ('Cc1nnc2csnc2n1',  '3-methylisothiazolo[3,4-e][1,2,4]triazine'),
    ('Cc1snc2ncnnc12',  '7-methylisothiazolo[3,4-e][1,2,4]triazine'),
    #
    # isothiazolo[4,3-e][1,2,4]triazine
    ('c1nnc2nscc2n1',  'isothiazolo[4,3-e][1,2,4]triazine'),
    ('Cc1nnc2nscc2n1',  '5-methylisothiazolo[4,3-e][1,2,4]triazine'),
    ('Cc1snc2nncnc12',  '3-methylisothiazolo[4,3-e][1,2,4]triazine'),
    #
    # isothiazolo[5,4-e][1,2,4]triazine
    ('c1nnc2cnsc2n1',  'isothiazolo[5,4-e][1,2,4]triazine'),
    ('Cc1nnc2cnsc2n1',  '3-methylisothiazolo[5,4-e][1,2,4]triazine'),
    ('Cc1nsc2ncnnc12',  '7-methylisothiazolo[5,4-e][1,2,4]triazine'),
    #
    # isoxazolo[3,4-e][1,2,4]triazine
    ('c1nnc2conc2n1',  'isoxazolo[3,4-e][1,2,4]triazine'),
    ('Cc1nnc2conc2n1',  '3-methylisoxazolo[3,4-e][1,2,4]triazine'),
    ('Cc1onc2ncnnc12',  '7-methylisoxazolo[3,4-e][1,2,4]triazine'),
    #
    # isoxazolo[4,3-e][1,2,4]triazine
    ('c1nnc2nocc2n1',  'isoxazolo[4,3-e][1,2,4]triazine'),
    ('Cc1nnc2nocc2n1',  '5-methylisoxazolo[4,3-e][1,2,4]triazine'),
    ('Cc1onc2nncnc12',  '3-methylisoxazolo[4,3-e][1,2,4]triazine'),
    #
    # isoxazolo[5,4-e][1,2,4]triazine
    ('c1nnc2cnoc2n1',  'isoxazolo[5,4-e][1,2,4]triazine'),
    ('Cc1nnc2cnoc2n1',  '3-methylisoxazolo[5,4-e][1,2,4]triazine'),
    ('Cc1noc2ncnnc12',  '7-methylisoxazolo[5,4-e][1,2,4]triazine'),
    #
    # isoxazolo[4,5-e][1,2,4]triazine
    ('c1nnc2oncc2n1',  'isoxazolo[4,5-e][1,2,4]triazine'),
    ('Cc1nnc2oncc2n1',  '5-methylisoxazolo[4,5-e][1,2,4]triazine'),
    ('Cc1noc2nncnc12',  '3-methylisoxazolo[4,5-e][1,2,4]triazine'),
    #
    # isothiazolo[4,5-e][1,2,4]triazine
    ('c1nnc2sncc2n1',  'isothiazolo[4,5-e][1,2,4]triazine'),
    ('Cc1nnc2sncc2n1',  '5-methylisothiazolo[4,5-e][1,2,4]triazine'),
    ('Cc1nsc2nncnc12',  '3-methylisothiazolo[4,5-e][1,2,4]triazine'),
    #
    # oxazolo[4,5-e][1,2,4]triazine
    ('c1nnc2ocnc2n1',  'oxazolo[4,5-e][1,2,4]triazine'),
    ('Cc1nnc2ocnc2n1',  '3-methyloxazolo[4,5-e][1,2,4]triazine'),
    ('Cc1nc2ncnnc2o1',  '6-methyloxazolo[4,5-e][1,2,4]triazine'),
    #
    # oxazolo[5,4-e][1,2,4]triazine
    ('c1nnc2ncoc2n1',  'oxazolo[5,4-e][1,2,4]triazine'),
    ('Cc1nnc2ncoc2n1',  '3-methyloxazolo[5,4-e][1,2,4]triazine'),
    ('Cc1nc2nncnc2o1',  '6-methyloxazolo[5,4-e][1,2,4]triazine'),
    #
    # thiazolo[4,5-e][1,2,4]triazine
    ('c1nnc2scnc2n1',  'thiazolo[4,5-e][1,2,4]triazine'),
    ('Cc1nnc2scnc2n1',  '3-methylthiazolo[4,5-e][1,2,4]triazine'),
    ('Cc1nc2ncnnc2s1',  '6-methylthiazolo[4,5-e][1,2,4]triazine'),
    #
    # thiazolo[5,4-e][1,2,4]triazine
    ('c1nnc2ncsc2n1',  'thiazolo[5,4-e][1,2,4]triazine'),
    ('Cc1nnc2ncsc2n1',  '3-methylthiazolo[5,4-e][1,2,4]triazine'),
    ('Cc1nc2nncnc2s1',  '6-methylthiazolo[5,4-e][1,2,4]triazine'),
])
def test_phase603(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
