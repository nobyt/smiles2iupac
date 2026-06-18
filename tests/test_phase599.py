"""Phase 599: isothiazolo fused with pyridine and pyridazine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isothiazolo[3,4-b]pyridine
    ('c1cnc2nscc2c1',  'isothiazolo[3,4-b]pyridine'),
    ('Cc1cnc2nscc2c1',  '5-methylisothiazolo[3,4-b]pyridine'),
    ('Cc1ccc2csnc2n1',  '6-methylisothiazolo[3,4-b]pyridine'),
    ('Cc1snc2ncccc12',  '3-methylisothiazolo[3,4-b]pyridine'),
    ('Cc1ccnc2nscc12',  '4-methylisothiazolo[3,4-b]pyridine'),
    #
    # isothiazolo[4,3-b]pyridine
    ('c1cnc2csnc2c1',  'isothiazolo[4,3-b]pyridine'),
    ('Cc1cnc2csnc2c1',  '6-methylisothiazolo[4,3-b]pyridine'),
    ('Cc1ccc2nscc2n1',  '5-methylisothiazolo[4,3-b]pyridine'),
    ('Cc1snc2cccnc12',  '3-methylisothiazolo[4,3-b]pyridine'),
    ('Cc1ccnc2csnc12',  '7-methylisothiazolo[4,3-b]pyridine'),
    #
    # isothiazolo[4,5-b]pyridine
    ('c1cnc2cnsc2c1',  'isothiazolo[4,5-b]pyridine'),
    ('Cc1cnc2cnsc2c1',  '6-methylisothiazolo[4,5-b]pyridine'),
    ('Cc1ccc2sncc2n1',  '5-methylisothiazolo[4,5-b]pyridine'),
    ('Cc1nsc2cccnc12',  '3-methylisothiazolo[4,5-b]pyridine'),
    ('Cc1ccnc2cnsc12',  '7-methylisothiazolo[4,5-b]pyridine'),
    #
    # isothiazolo[5,4-b]pyridine
    ('c1cnc2sncc2c1',  'isothiazolo[5,4-b]pyridine'),
    ('Cc1cnc2sncc2c1',  '5-methylisothiazolo[5,4-b]pyridine'),
    ('Cc1ccc2cnsc2n1',  '6-methylisothiazolo[5,4-b]pyridine'),
    ('Cc1nsc2ncccc12',  '3-methylisothiazolo[5,4-b]pyridine'),
    ('Cc1ccnc2sncc12',  '4-methylisothiazolo[5,4-b]pyridine'),
    #
    # isothiazolo[3,4-c]pyridine
    ('c1cc2csnc2cn1',  'isothiazolo[3,4-c]pyridine'),
    ('Cc1cc2csnc2cn1',  '5-methylisothiazolo[3,4-c]pyridine'),
    ('Cc1cncc2nscc12',  '4-methylisothiazolo[3,4-c]pyridine'),
    ('Cc1snc2cnccc12',  '3-methylisothiazolo[3,4-c]pyridine'),
    ('Cc1nccc2csnc12',  '7-methylisothiazolo[3,4-c]pyridine'),
    #
    # isothiazolo[4,3-c]pyridine
    ('c1cc2nscc2cn1',  'isothiazolo[4,3-c]pyridine'),
    ('Cc1cc2nscc2cn1',  '6-methylisothiazolo[4,3-c]pyridine'),
    ('Cc1cncc2csnc12',  '7-methylisothiazolo[4,3-c]pyridine'),
    ('Cc1snc2ccncc12',  '3-methylisothiazolo[4,3-c]pyridine'),
    ('Cc1nccc2nscc12',  '4-methylisothiazolo[4,3-c]pyridine'),
    #
    # isothiazolo[4,5-c]pyridine
    ('c1cc2sncc2cn1',  'isothiazolo[4,5-c]pyridine'),
    ('Cc1cc2sncc2cn1',  '6-methylisothiazolo[4,5-c]pyridine'),
    ('Cc1cncc2cnsc12',  '7-methylisothiazolo[4,5-c]pyridine'),
    ('Cc1nsc2ccncc12',  '3-methylisothiazolo[4,5-c]pyridine'),
    ('Cc1nccc2sncc12',  '4-methylisothiazolo[4,5-c]pyridine'),
    #
    # isothiazolo[5,4-c]pyridine
    ('c1cc2cnsc2cn1',  'isothiazolo[5,4-c]pyridine'),
    ('Cc1cc2cnsc2cn1',  '5-methylisothiazolo[5,4-c]pyridine'),
    ('Cc1cncc2sncc12',  '4-methylisothiazolo[5,4-c]pyridine'),
    ('Cc1nsc2cnccc12',  '3-methylisothiazolo[5,4-c]pyridine'),
    ('Cc1nccc2cnsc12',  '7-methylisothiazolo[5,4-c]pyridine'),
    #
    # isothiazolo[3,4-c]pyridazine
    ('c1cc2csnc2nn1',  'isothiazolo[3,4-c]pyridazine'),
    ('Cc1cc2csnc2nn1',  '5-methylisothiazolo[3,4-c]pyridazine'),
    ('Cc1cnnc2nscc12',  '4-methylisothiazolo[3,4-c]pyridazine'),
    ('Cc1snc2nnccc12',  '3-methylisothiazolo[3,4-c]pyridazine'),
    #
    # isothiazolo[4,3-c]pyridazine
    ('c1cc2nscc2nn1',  'isothiazolo[4,3-c]pyridazine'),
    ('Cc1cc2nscc2nn1',  '6-methylisothiazolo[4,3-c]pyridazine'),
    ('Cc1cnnc2csnc12',  '7-methylisothiazolo[4,3-c]pyridazine'),
    ('Cc1snc2ccnnc12',  '3-methylisothiazolo[4,3-c]pyridazine'),
    #
    # isothiazolo[4,5-c]pyridazine
    ('c1cc2sncc2nn1',  'isothiazolo[4,5-c]pyridazine'),
    ('Cc1cc2sncc2nn1',  '6-methylisothiazolo[4,5-c]pyridazine'),
    ('Cc1cnnc2cnsc12',  '7-methylisothiazolo[4,5-c]pyridazine'),
    ('Cc1nsc2ccnnc12',  '3-methylisothiazolo[4,5-c]pyridazine'),
    #
    # isothiazolo[5,4-c]pyridazine
    ('c1cc2cnsc2nn1',  'isothiazolo[5,4-c]pyridazine'),
    ('Cc1cc2cnsc2nn1',  '5-methylisothiazolo[5,4-c]pyridazine'),
    ('Cc1cnnc2sncc12',  '4-methylisothiazolo[5,4-c]pyridazine'),
    ('Cc1nsc2nnccc12',  '3-methylisothiazolo[5,4-c]pyridazine'),
])
def test_phase599(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
