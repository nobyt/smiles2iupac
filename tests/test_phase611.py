"""Phase 611: thieno/furo/isoxazolo fused with pyrimidine (d-fused)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[2,3-d]pyrimidine
    ('c1ncc2ccsc2n1',  'thieno[2,3-d]pyrimidine'),
    ('Cc1ncc2ccsc2n1',  '2-methylthieno[2,3-d]pyrimidine'),
    ('Cc1ncnc2sccc12',  '4-methylthieno[2,3-d]pyrimidine'),
    ('Cc1csc2ncncc12',  '5-methylthieno[2,3-d]pyrimidine'),
    ('Cc1cc2cncnc2s1',  '6-methylthieno[2,3-d]pyrimidine'),
    #
    # thieno[3,4-d]pyrimidine
    ('c1ncc2cscc2n1',  'thieno[3,4-d]pyrimidine'),
    ('Cc1ncc2cscc2n1',  '2-methylthieno[3,4-d]pyrimidine'),
    ('Cc1ncnc2cscc12',  '4-methylthieno[3,4-d]pyrimidine'),
    ('Cc1scc2ncncc12',  '5-methylthieno[3,4-d]pyrimidine'),
    ('Cc1scc2cncnc12',  '7-methylthieno[3,4-d]pyrimidine'),
    #
    # thieno[3,2-d]pyrimidine
    ('c1ncc2sccc2n1',  'thieno[3,2-d]pyrimidine'),
    ('Cc1ncc2sccc2n1',  '2-methylthieno[3,2-d]pyrimidine'),
    ('Cc1ncnc2ccsc12',  '4-methylthieno[3,2-d]pyrimidine'),
    ('Cc1cc2ncncc2s1',  '6-methylthieno[3,2-d]pyrimidine'),
    ('Cc1csc2cncnc12',  '7-methylthieno[3,2-d]pyrimidine'),
    #
    # furo[2,3-d]pyrimidine
    ('c1ncc2ccoc2n1',  'furo[2,3-d]pyrimidine'),
    ('Cc1ncc2ccoc2n1',  '2-methylfuro[2,3-d]pyrimidine'),
    ('Cc1ncnc2occc12',  '4-methylfuro[2,3-d]pyrimidine'),
    ('Cc1coc2ncncc12',  '5-methylfuro[2,3-d]pyrimidine'),
    ('Cc1cc2cncnc2o1',  '6-methylfuro[2,3-d]pyrimidine'),
    #
    # furo[3,4-d]pyrimidine
    ('c1ncc2cocc2n1',  'furo[3,4-d]pyrimidine'),
    ('Cc1ncc2cocc2n1',  '2-methylfuro[3,4-d]pyrimidine'),
    ('Cc1ncnc2cocc12',  '4-methylfuro[3,4-d]pyrimidine'),
    ('Cc1occ2ncncc12',  '5-methylfuro[3,4-d]pyrimidine'),
    ('Cc1occ2cncnc12',  '7-methylfuro[3,4-d]pyrimidine'),
    #
    # furo[3,2-d]pyrimidine
    ('c1ncc2occc2n1',  'furo[3,2-d]pyrimidine'),
    ('Cc1ncc2occc2n1',  '2-methylfuro[3,2-d]pyrimidine'),
    ('Cc1ncnc2ccoc12',  '4-methylfuro[3,2-d]pyrimidine'),
    ('Cc1cc2ncncc2o1',  '6-methylfuro[3,2-d]pyrimidine'),
    ('Cc1coc2cncnc12',  '7-methylfuro[3,2-d]pyrimidine'),
    #
    # isoxazolo[4,3-d]pyrimidine
    ('c1ncc2nocc2n1',  'isoxazolo[4,3-d]pyrimidine'),
    ('Cc1ncc2nocc2n1',  '5-methylisoxazolo[4,3-d]pyrimidine'),
    ('Cc1ncnc2conc12',  '7-methylisoxazolo[4,3-d]pyrimidine'),
    ('Cc1onc2cncnc12',  '3-methylisoxazolo[4,3-d]pyrimidine'),
    #
    # isoxazolo[4,5-d]pyrimidine
    ('c1ncc2oncc2n1',  'isoxazolo[4,5-d]pyrimidine'),
    ('Cc1ncc2oncc2n1',  '5-methylisoxazolo[4,5-d]pyrimidine'),
    ('Cc1ncnc2cnoc12',  '7-methylisoxazolo[4,5-d]pyrimidine'),
    ('Cc1noc2cncnc12',  '3-methylisoxazolo[4,5-d]pyrimidine'),
    #
    # isoxazolo[3,4-d]pyrimidine
    ('c1ncc2conc2n1',  'isoxazolo[3,4-d]pyrimidine'),
    ('Cc1ncc2conc2n1',  '6-methylisoxazolo[3,4-d]pyrimidine'),
    ('Cc1ncnc2nocc12',  '4-methylisoxazolo[3,4-d]pyrimidine'),
    ('Cc1onc2ncncc12',  '3-methylisoxazolo[3,4-d]pyrimidine'),
    #
    # isoxazolo[5,4-d]pyrimidine
    ('c1ncc2cnoc2n1',  'isoxazolo[5,4-d]pyrimidine'),
    ('Cc1ncc2cnoc2n1',  '6-methylisoxazolo[5,4-d]pyrimidine'),
    ('Cc1ncnc2oncc12',  '4-methylisoxazolo[5,4-d]pyrimidine'),
    ('Cc1noc2ncncc12',  '3-methylisoxazolo[5,4-d]pyrimidine'),
])
def test_phase611(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
