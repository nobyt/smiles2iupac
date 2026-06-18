"""Phase 610: oxazolo/thiazolo-d-pyrimidine; isoxazolo-b-pyridine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # oxazolo[5,4-d]pyrimidine
    ('c1ncc2ncoc2n1',  'oxazolo[5,4-d]pyrimidine'),
    ('Cc1ncc2ncoc2n1',  '5-methyloxazolo[5,4-d]pyrimidine'),
    ('Cc1ncnc2ocnc12',  '7-methyloxazolo[5,4-d]pyrimidine'),
    ('Cc1nc2cncnc2o1',  '2-methyloxazolo[5,4-d]pyrimidine'),
    #
    # oxazolo[4,5-d]pyrimidine
    ('c1ncc2ocnc2n1',  'oxazolo[4,5-d]pyrimidine'),
    ('Cc1ncc2ocnc2n1',  '5-methyloxazolo[4,5-d]pyrimidine'),
    ('Cc1ncnc2ncoc12',  '7-methyloxazolo[4,5-d]pyrimidine'),
    ('Cc1nc2ncncc2o1',  '2-methyloxazolo[4,5-d]pyrimidine'),
    #
    # thiazolo[5,4-d]pyrimidine
    ('c1ncc2ncsc2n1',  'thiazolo[5,4-d]pyrimidine'),
    ('Cc1ncc2ncsc2n1',  '5-methylthiazolo[5,4-d]pyrimidine'),
    ('Cc1ncnc2scnc12',  '7-methylthiazolo[5,4-d]pyrimidine'),
    ('Cc1nc2cncnc2s1',  '2-methylthiazolo[5,4-d]pyrimidine'),
    #
    # thiazolo[4,5-d]pyrimidine
    ('c1ncc2scnc2n1',  'thiazolo[4,5-d]pyrimidine'),
    ('Cc1ncc2scnc2n1',  '5-methylthiazolo[4,5-d]pyrimidine'),
    ('Cc1ncnc2ncsc12',  '7-methylthiazolo[4,5-d]pyrimidine'),
    ('Cc1nc2ncncc2s1',  '2-methylthiazolo[4,5-d]pyrimidine'),
    #
    # isoxazolo[3,4-b]pyridine
    ('c1cnc2nocc2c1',  'isoxazolo[3,4-b]pyridine'),
    ('Cc1cnc2nocc2c1',  '5-methylisoxazolo[3,4-b]pyridine'),
    ('Cc1ccc2conc2n1',  '6-methylisoxazolo[3,4-b]pyridine'),
    ('Cc1onc2ncccc12',  '3-methylisoxazolo[3,4-b]pyridine'),
    ('Cc1ccnc2nocc12',  '4-methylisoxazolo[3,4-b]pyridine'),
    #
    # isoxazolo[4,5-b]pyridine
    ('c1cnc2cnoc2c1',  'isoxazolo[4,5-b]pyridine'),
    ('Cc1cnc2cnoc2c1',  '6-methylisoxazolo[4,5-b]pyridine'),
    ('Cc1ccc2oncc2n1',  '5-methylisoxazolo[4,5-b]pyridine'),
    ('Cc1noc2cccnc12',  '3-methylisoxazolo[4,5-b]pyridine'),
    ('Cc1ccnc2cnoc12',  '7-methylisoxazolo[4,5-b]pyridine'),
    #
    # isoxazolo[5,4-b]pyridine
    ('c1cnc2oncc2c1',  'isoxazolo[5,4-b]pyridine'),
    ('Cc1cnc2oncc2c1',  '5-methylisoxazolo[5,4-b]pyridine'),
    ('Cc1ccc2cnoc2n1',  '6-methylisoxazolo[5,4-b]pyridine'),
    ('Cc1noc2ncccc12',  '3-methylisoxazolo[5,4-b]pyridine'),
    ('Cc1ccnc2oncc12',  '4-methylisoxazolo[5,4-b]pyridine'),
    #
    # isoxazolo[4,3-b]pyridine
    ('c1cnc2conc2c1',  'isoxazolo[4,3-b]pyridine'),
    ('Cc1cnc2conc2c1',  '6-methylisoxazolo[4,3-b]pyridine'),
    ('Cc1ccc2nocc2n1',  '5-methylisoxazolo[4,3-b]pyridine'),
    ('Cc1onc2cccnc12',  '3-methylisoxazolo[4,3-b]pyridine'),
    ('Cc1ccnc2conc12',  '7-methylisoxazolo[4,3-b]pyridine'),
])
def test_phase610(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
