"""Phase 600: oxazolo/thiazolo/thieno fused with pyridine/pyrazine; imidazo[1,2-c]pyrimidine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thieno[3,4-b]pyridine
    ('c1cnc2cscc2c1',  'thieno[3,4-b]pyridine'),
    ('Cc1cnc2cscc2c1',  '3-methylthieno[3,4-b]pyridine'),
    ('Cc1ccc2cscc2n1',  '2-methylthieno[3,4-b]pyridine'),
    ('Cc1scc2cccnc12',  '7-methylthieno[3,4-b]pyridine'),
    ('Cc1scc2ncccc12',  '5-methylthieno[3,4-b]pyridine'),
    ('Cc1ccnc2cscc12',  '4-methylthieno[3,4-b]pyridine'),
    #
    # oxazolo[5,4-b]pyridine
    ('c1cnc2ocnc2c1',  'oxazolo[5,4-b]pyridine'),
    ('Cc1cnc2ocnc2c1',  '6-methyloxazolo[5,4-b]pyridine'),
    ('Cc1ccc2ncoc2n1',  '5-methyloxazolo[5,4-b]pyridine'),
    ('Cc1nc2cccnc2o1',  '2-methyloxazolo[5,4-b]pyridine'),
    ('Cc1ccnc2ocnc12',  '7-methyloxazolo[5,4-b]pyridine'),
    #
    # oxazolo[4,5-b]pyridine
    ('c1cnc2ncoc2c1',  'oxazolo[4,5-b]pyridine'),
    ('Cc1cnc2ncoc2c1',  '6-methyloxazolo[4,5-b]pyridine'),
    ('Cc1ccc2ocnc2n1',  '5-methyloxazolo[4,5-b]pyridine'),
    ('Cc1nc2ncccc2o1',  '2-methyloxazolo[4,5-b]pyridine'),
    ('Cc1ccnc2ncoc12',  '7-methyloxazolo[4,5-b]pyridine'),
    #
    # thiazolo[5,4-b]pyridine
    ('c1cnc2scnc2c1',  'thiazolo[5,4-b]pyridine'),
    ('Cc1cnc2scnc2c1',  '6-methylthiazolo[5,4-b]pyridine'),
    ('Cc1ccc2ncsc2n1',  '5-methylthiazolo[5,4-b]pyridine'),
    ('Cc1nc2cccnc2s1',  '2-methylthiazolo[5,4-b]pyridine'),
    ('Cc1ccnc2scnc12',  '7-methylthiazolo[5,4-b]pyridine'),
    #
    # thiazolo[4,5-b]pyridine
    ('c1cnc2ncsc2c1',  'thiazolo[4,5-b]pyridine'),
    ('Cc1cnc2ncsc2c1',  '6-methylthiazolo[4,5-b]pyridine'),
    ('Cc1ccc2scnc2n1',  '5-methylthiazolo[4,5-b]pyridine'),
    ('Cc1nc2ncccc2s1',  '2-methylthiazolo[4,5-b]pyridine'),
    ('Cc1ccnc2ncsc12',  '7-methylthiazolo[4,5-b]pyridine'),
    #
    # oxazolo[4,5-e]pyrazine
    ('c1cnc2ocnc2n1',  'oxazolo[4,5-e]pyrazine'),
    ('Cc1cnc2ocnc2n1',  '5-methyloxazolo[4,5-e]pyrazine'),
    ('Cc1cnc2ncoc2n1',  '6-methyloxazolo[4,5-e]pyrazine'),
    ('Cc1nc2nccnc2o1',  '2-methyloxazolo[4,5-e]pyrazine'),
    #
    # thiazolo[4,5-e]pyrazine
    ('c1cnc2scnc2n1',  'thiazolo[4,5-e]pyrazine'),
    ('Cc1cnc2scnc2n1',  '5-methylthiazolo[4,5-e]pyrazine'),
    ('Cc1cnc2ncsc2n1',  '6-methylthiazolo[4,5-e]pyrazine'),
    ('Cc1nc2nccnc2s1',  '2-methylthiazolo[4,5-e]pyrazine'),
    #
    # thieno[3,4-c]pyridine
    ('c1cc2cscc2cn1',  'thieno[3,4-c]pyridine'),
    ('Cc1cc2cscc2cn1',  '6-methylthieno[3,4-c]pyridine'),
    ('Cc1cncc2cscc12',  '7-methylthieno[3,4-c]pyridine'),
    ('Cc1scc2cnccc12',  '1-methylthieno[3,4-c]pyridine'),
    ('Cc1scc2ccncc12',  '3-methylthieno[3,4-c]pyridine'),
    ('Cc1nccc2cscc12',  '4-methylthieno[3,4-c]pyridine'),
    #
    # imidazo[1,2-c]pyrimidine
    ('c1cc2nccn2cn1',  'imidazo[1,2-c]pyrimidine'),
    ('Cc1cc2nccn2cn1',  '7-methylimidazo[1,2-c]pyrimidine'),
    ('Cc1cncn2ccnc12',  '8-methylimidazo[1,2-c]pyrimidine'),
    ('Cc1cn2cnccc2n1',  '2-methylimidazo[1,2-c]pyrimidine'),
    ('Cc1cnc2ccncn12',  '3-methylimidazo[1,2-c]pyrimidine'),
    ('Cc1nccc2nccn12',  '5-methylimidazo[1,2-c]pyrimidine'),
])
def test_phase600(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
