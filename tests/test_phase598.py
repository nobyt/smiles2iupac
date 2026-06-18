"""Phase 598: thieno/furo-pyrazine, [1,2,3]triazolo-pyrimidine, isoxazolo/isothiazolo-pyrimidine/pyrazine (5+6 bicyclics)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[2,3-e]pyrazine
    ('c1cnc2sccc2n1',  'thieno[2,3-e]pyrazine'),
    ('Cc1cnc2sccc2n1',  '2-methylthieno[2,3-e]pyrazine'),
    ('Cc1cnc2ccsc2n1',  '3-methylthieno[2,3-e]pyrazine'),
    ('Cc1cc2nccnc2s1',  '6-methylthieno[2,3-e]pyrazine'),
    ('Cc1csc2nccnc12',  '7-methylthieno[2,3-e]pyrazine'),
    #
    # furo[2,3-e]pyrazine
    ('c1cnc2occc2n1',  'furo[2,3-e]pyrazine'),
    ('Cc1cnc2occc2n1',  '2-methylfuro[2,3-e]pyrazine'),
    ('Cc1cnc2ccoc2n1',  '3-methylfuro[2,3-e]pyrazine'),
    ('Cc1cc2nccnc2o1',  '6-methylfuro[2,3-e]pyrazine'),
    ('Cc1coc2nccnc12',  '7-methylfuro[2,3-e]pyrazine'),
    #
    # furo[3,4-e]pyrazine
    ('c1cnc2cocc2n1',  'furo[3,4-e]pyrazine'),
    ('Cc1cnc2cocc2n1',  '2-methylfuro[3,4-e]pyrazine'),
    ('Cc1occ2nccnc12',  '5-methylfuro[3,4-e]pyrazine'),
    #
    # furo[2,3-b]pyridine
    ('c1cnc2occc2c1',  'furo[2,3-b]pyridine'),
    ('Cc1cnc2occc2c1',  '5-methylfuro[2,3-b]pyridine'),
    ('Cc1ccc2ccoc2n1',  '6-methylfuro[2,3-b]pyridine'),
    ('Cc1cc2cccnc2o1',  '2-methylfuro[2,3-b]pyridine'),
    ('Cc1coc2ncccc12',  '3-methylfuro[2,3-b]pyridine'),
    ('Cc1ccnc2occc12',  '4-methylfuro[2,3-b]pyridine'),
    #
    # thieno[3,4-e]pyrazine
    ('c1cnc2cscc2n1',  'thieno[3,4-e]pyrazine'),
    ('Cc1cnc2cscc2n1',  '2-methylthieno[3,4-e]pyrazine'),
    ('Cc1scc2nccnc12',  '5-methylthieno[3,4-e]pyrazine'),
    #
    # [1,2,3]triazolo[1,5-a]pyrimidine
    ('c1cnc2cnnn2c1',  '[1,2,3]triazolo[1,5-a]pyrimidine'),
    ('Cc1cnc2cnnn2c1',  '6-methyl[1,2,3]triazolo[1,5-a]pyrimidine'),
    ('Cc1ccn2nncc2n1',  '5-methyl[1,2,3]triazolo[1,5-a]pyrimidine'),
    ('Cc1nnn2cccnc12',  '3-methyl[1,2,3]triazolo[1,5-a]pyrimidine'),
    ('Cc1ccnc2cnnn12',  '7-methyl[1,2,3]triazolo[1,5-a]pyrimidine'),
    #
    # isoxazolo[3,4-e]pyrazine
    ('c1cnc2nocc2n1',  'isoxazolo[3,4-e]pyrazine'),
    ('Cc1cnc2nocc2n1',  '5-methylisoxazolo[3,4-e]pyrazine'),
    ('Cc1cnc2conc2n1',  '6-methylisoxazolo[3,4-e]pyrazine'),
    ('Cc1onc2nccnc12',  '3-methylisoxazolo[3,4-e]pyrazine'),
    #
    # isoxazolo[4,5-e]pyrazine
    ('c1cnc2oncc2n1',  'isoxazolo[4,5-e]pyrazine'),
    ('Cc1cnc2oncc2n1',  '5-methylisoxazolo[4,5-e]pyrazine'),
    ('Cc1cnc2cnoc2n1',  '6-methylisoxazolo[4,5-e]pyrazine'),
    ('Cc1noc2nccnc12',  '3-methylisoxazolo[4,5-e]pyrazine'),
    #
    # isothiazolo[4,5-d]pyrimidine
    ('c1ncc2sncc2n1',  'isothiazolo[4,5-d]pyrimidine'),
    ('Cc1ncc2sncc2n1',  '5-methylisothiazolo[4,5-d]pyrimidine'),
    ('Cc1ncnc2cnsc12',  '7-methylisothiazolo[4,5-d]pyrimidine'),
    ('Cc1nsc2cncnc12',  '3-methylisothiazolo[4,5-d]pyrimidine'),
    #
    # isothiazolo[4,3-d]pyrimidine
    ('c1ncc2nscc2n1',  'isothiazolo[4,3-d]pyrimidine'),
    ('Cc1ncc2nscc2n1',  '5-methylisothiazolo[4,3-d]pyrimidine'),
    ('Cc1ncnc2csnc12',  '7-methylisothiazolo[4,3-d]pyrimidine'),
    ('Cc1snc2cncnc12',  '3-methylisothiazolo[4,3-d]pyrimidine'),
    #
    # isothiazolo[5,4-d]pyrimidine
    ('c1ncc2cnsc2n1',  'isothiazolo[5,4-d]pyrimidine'),
    ('Cc1ncc2cnsc2n1',  '6-methylisothiazolo[5,4-d]pyrimidine'),
    ('Cc1ncnc2sncc12',  '4-methylisothiazolo[5,4-d]pyrimidine'),
    ('Cc1nsc2ncncc12',  '3-methylisothiazolo[5,4-d]pyrimidine'),
    #
    # isothiazolo[3,4-d]pyrimidine
    ('c1ncc2csnc2n1',  'isothiazolo[3,4-d]pyrimidine'),
    ('Cc1ncc2csnc2n1',  '6-methylisothiazolo[3,4-d]pyrimidine'),
    ('Cc1ncnc2nscc12',  '4-methylisothiazolo[3,4-d]pyrimidine'),
    ('Cc1snc2ncncc12',  '3-methylisothiazolo[3,4-d]pyrimidine'),
    #
    # isothiazolo[4,5-e]pyrazine
    ('c1cnc2sncc2n1',  'isothiazolo[4,5-e]pyrazine'),
    ('Cc1cnc2sncc2n1',  '5-methylisothiazolo[4,5-e]pyrazine'),
    ('Cc1cnc2cnsc2n1',  '6-methylisothiazolo[4,5-e]pyrazine'),
    ('Cc1nsc2nccnc12',  '3-methylisothiazolo[4,5-e]pyrazine'),
    #
    # isothiazolo[3,4-e]pyrazine
    ('c1cnc2nscc2n1',  'isothiazolo[3,4-e]pyrazine'),
    ('Cc1cnc2nscc2n1',  '5-methylisothiazolo[3,4-e]pyrazine'),
    ('Cc1cnc2csnc2n1',  '6-methylisothiazolo[3,4-e]pyrazine'),
    ('Cc1snc2nccnc12',  '3-methylisothiazolo[3,4-e]pyrazine'),
])
def test_phase598(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
