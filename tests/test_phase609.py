"""Phase 609: thiazolo/oxazolo/isoxazolo fused with pyridine and pyridazine (c-fused)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thiazolo[5,4-c]pyridine
    ('c1cc2ncsc2cn1',  'thiazolo[5,4-c]pyridine'),
    ('Cc1cc2ncsc2cn1',  '6-methylthiazolo[5,4-c]pyridine'),
    ('Cc1cncc2scnc12',  '7-methylthiazolo[5,4-c]pyridine'),
    ('Cc1nc2ccncc2s1',  '2-methylthiazolo[5,4-c]pyridine'),
    ('Cc1nccc2ncsc12',  '4-methylthiazolo[5,4-c]pyridine'),
    #
    # oxazolo[5,4-c]pyridine
    ('c1cc2ncoc2cn1',  'oxazolo[5,4-c]pyridine'),
    ('Cc1cc2ncoc2cn1',  '6-methyloxazolo[5,4-c]pyridine'),
    ('Cc1cncc2ocnc12',  '7-methyloxazolo[5,4-c]pyridine'),
    ('Cc1nc2ccncc2o1',  '2-methyloxazolo[5,4-c]pyridine'),
    ('Cc1nccc2ncoc12',  '4-methyloxazolo[5,4-c]pyridine'),
    #
    # thiazolo[4,5-c]pyridine
    ('c1cc2scnc2cn1',  'thiazolo[4,5-c]pyridine'),
    ('Cc1cc2scnc2cn1',  '6-methylthiazolo[4,5-c]pyridine'),
    ('Cc1cncc2ncsc12',  '7-methylthiazolo[4,5-c]pyridine'),
    ('Cc1nc2cnccc2s1',  '2-methylthiazolo[4,5-c]pyridine'),
    ('Cc1nccc2scnc12',  '4-methylthiazolo[4,5-c]pyridine'),
    #
    # thiazolo[4,5-c]pyridazine
    ('c1cc2scnc2nn1',  'thiazolo[4,5-c]pyridazine'),
    ('Cc1cc2scnc2nn1',  '3-methylthiazolo[4,5-c]pyridazine'),
    ('Cc1cnnc2ncsc12',  '4-methylthiazolo[4,5-c]pyridazine'),
    ('Cc1nc2nnccc2s1',  '6-methylthiazolo[4,5-c]pyridazine'),
    #
    # thiazolo[5,4-c]pyridazine
    ('c1cc2ncsc2nn1',  'thiazolo[5,4-c]pyridazine'),
    ('Cc1cc2ncsc2nn1',  '3-methylthiazolo[5,4-c]pyridazine'),
    ('Cc1cnnc2scnc12',  '4-methylthiazolo[5,4-c]pyridazine'),
    ('Cc1nc2ccnnc2s1',  '6-methylthiazolo[5,4-c]pyridazine'),
    #
    # oxazolo[4,5-c]pyridine
    ('c1cc2ocnc2cn1',  'oxazolo[4,5-c]pyridine'),
    ('Cc1cc2ocnc2cn1',  '6-methyloxazolo[4,5-c]pyridine'),
    ('Cc1cncc2ncoc12',  '7-methyloxazolo[4,5-c]pyridine'),
    ('Cc1nc2cnccc2o1',  '2-methyloxazolo[4,5-c]pyridine'),
    ('Cc1nccc2ocnc12',  '4-methyloxazolo[4,5-c]pyridine'),
    #
    # isoxazolo[5,4-c]pyridine
    ('c1cc2cnoc2cn1',  'isoxazolo[5,4-c]pyridine'),
    ('Cc1cc2cnoc2cn1',  '5-methylisoxazolo[5,4-c]pyridine'),
    ('Cc1cncc2oncc12',  '4-methylisoxazolo[5,4-c]pyridine'),
    ('Cc1noc2cnccc12',  '3-methylisoxazolo[5,4-c]pyridine'),
    ('Cc1nccc2cnoc12',  '7-methylisoxazolo[5,4-c]pyridine'),
    #
    # isoxazolo[4,5-c]pyridine
    ('c1cc2oncc2cn1',  'isoxazolo[4,5-c]pyridine'),
    ('Cc1cc2oncc2cn1',  '6-methylisoxazolo[4,5-c]pyridine'),
    ('Cc1cncc2cnoc12',  '7-methylisoxazolo[4,5-c]pyridine'),
    ('Cc1noc2ccncc12',  '3-methylisoxazolo[4,5-c]pyridine'),
    ('Cc1nccc2oncc12',  '4-methylisoxazolo[4,5-c]pyridine'),
    #
    # isoxazolo[4,3-c]pyridine
    ('c1cc2nocc2cn1',  'isoxazolo[4,3-c]pyridine'),
    ('Cc1cc2nocc2cn1',  '6-methylisoxazolo[4,3-c]pyridine'),
    ('Cc1cncc2conc12',  '7-methylisoxazolo[4,3-c]pyridine'),
    ('Cc1onc2ccncc12',  '3-methylisoxazolo[4,3-c]pyridine'),
    ('Cc1nccc2nocc12',  '4-methylisoxazolo[4,3-c]pyridine'),
    #
    # oxazolo[5,4-c]pyridazine
    ('c1cc2ncoc2nn1',  'oxazolo[5,4-c]pyridazine'),
    ('Cc1cc2ncoc2nn1',  '3-methyloxazolo[5,4-c]pyridazine'),
    ('Cc1cnnc2ocnc12',  '4-methyloxazolo[5,4-c]pyridazine'),
    ('Cc1nc2ccnnc2o1',  '6-methyloxazolo[5,4-c]pyridazine'),
    #
    # oxazolo[4,5-c]pyridazine
    ('c1cc2ocnc2nn1',  'oxazolo[4,5-c]pyridazine'),
    ('Cc1cc2ocnc2nn1',  '3-methyloxazolo[4,5-c]pyridazine'),
    ('Cc1cnnc2ncoc12',  '4-methyloxazolo[4,5-c]pyridazine'),
    ('Cc1nc2nnccc2o1',  '6-methyloxazolo[4,5-c]pyridazine'),
    #
    # isoxazolo[5,4-c]pyridazine
    ('c1cc2cnoc2nn1',  'isoxazolo[5,4-c]pyridazine'),
    ('Cc1cc2cnoc2nn1',  '5-methylisoxazolo[5,4-c]pyridazine'),
    ('Cc1cnnc2oncc12',  '4-methylisoxazolo[5,4-c]pyridazine'),
    ('Cc1noc2nnccc12',  '3-methylisoxazolo[5,4-c]pyridazine'),
    #
    # isoxazolo[4,5-c]pyridazine
    ('c1cc2oncc2nn1',  'isoxazolo[4,5-c]pyridazine'),
    ('Cc1cc2oncc2nn1',  '6-methylisoxazolo[4,5-c]pyridazine'),
    ('Cc1cnnc2cnoc12',  '7-methylisoxazolo[4,5-c]pyridazine'),
    ('Cc1noc2ccnnc12',  '3-methylisoxazolo[4,5-c]pyridazine'),
    #
    # isoxazolo[4,3-c]pyridazine
    ('c1cc2nocc2nn1',  'isoxazolo[4,3-c]pyridazine'),
    ('Cc1cc2nocc2nn1',  '6-methylisoxazolo[4,3-c]pyridazine'),
    ('Cc1cnnc2conc12',  '7-methylisoxazolo[4,3-c]pyridazine'),
    ('Cc1onc2ccnnc12',  '3-methylisoxazolo[4,3-c]pyridazine'),
    #
    # isoxazolo[3,4-c]pyridine
    ('c1cc2conc2cn1',  'isoxazolo[3,4-c]pyridine'),
    ('Cc1cc2conc2cn1',  '5-methylisoxazolo[3,4-c]pyridine'),
    ('Cc1cncc2nocc12',  '4-methylisoxazolo[3,4-c]pyridine'),
    ('Cc1onc2cnccc12',  '3-methylisoxazolo[3,4-c]pyridine'),
    ('Cc1nccc2conc12',  '7-methylisoxazolo[3,4-c]pyridine'),
    #
    # isoxazolo[3,4-c]pyridazine
    ('c1cc2conc2nn1',  'isoxazolo[3,4-c]pyridazine'),
    ('Cc1cc2conc2nn1',  '5-methylisoxazolo[3,4-c]pyridazine'),
    ('Cc1cnnc2nocc12',  '4-methylisoxazolo[3,4-c]pyridazine'),
    ('Cc1onc2nnccc12',  '3-methylisoxazolo[3,4-c]pyridazine'),
])
def test_phase609(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
