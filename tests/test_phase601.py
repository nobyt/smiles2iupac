"""Phase 601: pyrazolo/triazolo/imidazo fused with pyrazine/[1,2,4]triazine (N-bridged)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrazolo[1,5-a]pyrazine
    ('c1cn2nccc2cn1',  'pyrazolo[1,5-a]pyrazine'),
    ('Cc1cn2nccc2cn1',  '6-methylpyrazolo[1,5-a]pyrazine'),
    ('Cc1cncc2ccnn12',  '7-methylpyrazolo[1,5-a]pyrazine'),
    ('Cc1cc2cnccn2n1',  '2-methylpyrazolo[1,5-a]pyrazine'),
    ('Cc1cnn2ccncc12',  '3-methylpyrazolo[1,5-a]pyrazine'),
    ('Cc1nccn2nccc12',  '4-methylpyrazolo[1,5-a]pyrazine'),
    #
    # [1,2,4]triazolo[4,3-b][1,2,4]triazine
    ('c1cnn2cnnc2n1',  '[1,2,4]triazolo[4,3-b][1,2,4]triazine'),
    ('Cc1cnn2cnnc2n1',  '7-methyl[1,2,4]triazolo[4,3-b][1,2,4]triazine'),
    ('Cc1cnc2nncn2n1',  '6-methyl[1,2,4]triazolo[4,3-b][1,2,4]triazine'),
    ('Cc1nnc2nccnn12',  '3-methyl[1,2,4]triazolo[4,3-b][1,2,4]triazine'),
    #
    # [1,2,4]triazolo[1,5-a]pyrazine
    ('c1cn2ncnc2cn1',  '[1,2,4]triazolo[1,5-a]pyrazine'),
    ('Cc1cn2ncnc2cn1',  '6-methyl[1,2,4]triazolo[1,5-a]pyrazine'),
    ('Cc1cncc2ncnn12',  '5-methyl[1,2,4]triazolo[1,5-a]pyrazine'),
    ('Cc1nc2cnccn2n1',  '2-methyl[1,2,4]triazolo[1,5-a]pyrazine'),
    ('Cc1nccn2ncnc12',  '8-methyl[1,2,4]triazolo[1,5-a]pyrazine'),
    #
    # imidazo[1,5-b][1,2,4]triazine
    ('c1cnn2cncc2n1',  'imidazo[1,5-b][1,2,4]triazine'),
    ('Cc1cnn2cncc2n1',  '2-methylimidazo[1,5-b][1,2,4]triazine'),
    ('Cc1cnc2cncn2n1',  '3-methylimidazo[1,5-b][1,2,4]triazine'),
    ('Cc1ncc2nccnn12',  '6-methylimidazo[1,5-b][1,2,4]triazine'),
    ('Cc1ncn2nccnc12',  '8-methylimidazo[1,5-b][1,2,4]triazine'),
    #
    # imidazo[3,2-b][1,2,4]triazine
    ('c1cnn2ccnc2n1',  'imidazo[3,2-b][1,2,4]triazine'),
    ('Cc1cnn2ccnc2n1',  '3-methylimidazo[3,2-b][1,2,4]triazine'),
    ('Cc1cnc2nccn2n1',  '2-methylimidazo[3,2-b][1,2,4]triazine'),
    ('Cc1cnc2nccnn12',  '7-methylimidazo[3,2-b][1,2,4]triazine'),
    ('Cc1cn2nccnc2n1',  '6-methylimidazo[3,2-b][1,2,4]triazine'),
    #
    # imidazo[1,5-a]pyrazine
    ('c1cn2cncc2cn1',  'imidazo[1,5-a]pyrazine'),
    ('Cc1cn2cncc2cn1',  '6-methylimidazo[1,5-a]pyrazine'),
    ('Cc1cncc2cncn12',  '5-methylimidazo[1,5-a]pyrazine'),
    ('Cc1ncc2cnccn12',  '3-methylimidazo[1,5-a]pyrazine'),
    ('Cc1ncn2ccncc12',  '1-methylimidazo[1,5-a]pyrazine'),
    ('Cc1nccn2cncc12',  '8-methylimidazo[1,5-a]pyrazine'),
    #
    # [1,2,4]triazolo[4,3-a]pyrazine
    ('c1cn2cnnc2cn1',  '[1,2,4]triazolo[4,3-a]pyrazine'),
    ('Cc1cn2cnnc2cn1',  '6-methyl[1,2,4]triazolo[4,3-a]pyrazine'),
    ('Cc1cncc2nncn12',  '5-methyl[1,2,4]triazolo[4,3-a]pyrazine'),
    ('Cc1nnc2cnccn12',  '3-methyl[1,2,4]triazolo[4,3-a]pyrazine'),
    ('Cc1nccn2cnnc12',  '8-methyl[1,2,4]triazolo[4,3-a]pyrazine'),
    #
    # pyrazolo[1,5-b][1,2,4]triazine
    ('c1cnn2nccc2n1',  'pyrazolo[1,5-b][1,2,4]triazine'),
    ('Cc1cnn2nccc2n1',  '2-methylpyrazolo[1,5-b][1,2,4]triazine'),
    ('Cc1cnc2ccnn2n1',  '3-methylpyrazolo[1,5-b][1,2,4]triazine'),
    ('Cc1cc2nccnn2n1',  '7-methylpyrazolo[1,5-b][1,2,4]triazine'),
    ('Cc1cnn2nccnc12',  '8-methylpyrazolo[1,5-b][1,2,4]triazine'),
])
def test_phase601(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
