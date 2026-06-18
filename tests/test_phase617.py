"""Phase 617: NH-containing bicyclic heteroaromatics — methyl-substituted pyrrolo, pyrazolo,
imidazo, triazolo fused systems and purines."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-pyrrolo[2,3-b]pyridine
    ('c1cnc2[nH]ccc2c1',    '1H-pyrrolo[2,3-b]pyridine'),
    ('Cc1cc2cccnc2[nH]1',   '2-methyl-1H-pyrrolo[2,3-b]pyridine'),
    ('Cc1c[nH]c2ncccc12',   '3-methyl-1H-pyrrolo[2,3-b]pyridine'),
    ('Cc1ccnc2[nH]ccc12',   '4-methyl-1H-pyrrolo[2,3-b]pyridine'),
    ('Cc1cnc2[nH]ccc2c1',   '5-methyl-1H-pyrrolo[2,3-b]pyridine'),
    ('Cc1ccc2cc[nH]c2n1',   '6-methyl-1H-pyrrolo[2,3-b]pyridine'),
    # 1H-pyrrolo[3,2-b]pyridine
    ('c1cnc2cc[nH]c2c1',    '1H-pyrrolo[3,2-b]pyridine'),
    ('Cc1cc2ncccc2[nH]1',   '2-methyl-1H-pyrrolo[3,2-b]pyridine'),
    ('Cc1c[nH]c2cccnc12',   '3-methyl-1H-pyrrolo[3,2-b]pyridine'),
    ('Cc1ccc2[nH]ccc2n1',   '5-methyl-1H-pyrrolo[3,2-b]pyridine'),
    ('Cc1cnc2cc[nH]c2c1',   '6-methyl-1H-pyrrolo[3,2-b]pyridine'),
    ('Cc1ccnc2cc[nH]c12',   '7-methyl-1H-pyrrolo[3,2-b]pyridine'),
    # 1H-pyrazolo[4,5-b]pyridine
    ('c1cnc2cn[nH]c2c1',    '1H-pyrazolo[4,5-b]pyridine'),
    ('Cc1n[nH]c2cccnc12',   '3-methyl-1H-pyrazolo[4,5-b]pyridine'),
    ('Cc1ccc2[nH]ncc2n1',   '5-methyl-1H-pyrazolo[4,5-b]pyridine'),
    ('Cc1cnc2cn[nH]c2c1',   '6-methyl-1H-pyrazolo[4,5-b]pyridine'),
    ('Cc1ccnc2cn[nH]c12',   '7-methyl-1H-pyrazolo[4,5-b]pyridine'),
    # 1H-pyrazolo[3,4-b]pyridine
    ('c1cnc2[nH]ncc2c1',    '1H-pyrazolo[3,4-b]pyridine'),
    ('Cc1n[nH]c2ncccc12',   '3-methyl-1H-pyrazolo[3,4-b]pyridine'),
    ('Cc1ccnc2[nH]ncc12',   '4-methyl-1H-pyrazolo[3,4-b]pyridine'),
    ('Cc1cnc2[nH]ncc2c1',   '5-methyl-1H-pyrazolo[3,4-b]pyridine'),
    ('Cc1ccc2cn[nH]c2n1',   '6-methyl-1H-pyrazolo[3,4-b]pyridine'),
    # 1H-imidazo[4,5-b]pyridine
    ('c1cnc2[nH]cnc2c1',    '1H-imidazo[4,5-b]pyridine'),
    ('Cc1nc2cccnc2[nH]1',   '2-methyl-1H-imidazo[4,5-b]pyridine'),
    ('Cc1ccc2nc[nH]c2n1',   '5-methyl-1H-imidazo[4,5-b]pyridine'),
    ('Cc1cnc2[nH]cnc2c1',   '6-methyl-1H-imidazo[4,5-b]pyridine'),
    ('Cc1ccnc2[nH]cnc12',   '7-methyl-1H-imidazo[4,5-b]pyridine'),
    # 7H-purine
    ('c1ncc2[nH]cnc2n1',    '7H-purine'),
    ('Cc1ncc2[nH]cnc2n1',   '2-methyl-7H-purine'),
    ('Cc1ncnc2nc[nH]c12',   '6-methyl-7H-purine'),
    ('Cc1nc2ncncc2[nH]1',   '8-methyl-7H-purine'),
    # 9H-purine
    ('c1ncc2nc[nH]c2n1',    '9H-purine'),
    ('Cc1ncc2nc[nH]c2n1',   '2-methyl-9H-purine'),
    ('Cc1ncnc2[nH]cnc12',   '6-methyl-9H-purine'),
    ('Cc1nc2cncnc2[nH]1',   '8-methyl-9H-purine'),
    # 2H-[1,2,3]triazolo[4,5-b]pyridine
    ('c1cnc2n[nH]nc2c1',    '2H-[1,2,3]triazolo[4,5-b]pyridine'),
    ('Cc1ccc2n[nH]nc2n1',   '5-methyl-2H-[1,2,3]triazolo[4,5-b]pyridine'),
    ('Cc1cnc2n[nH]nc2c1',   '6-methyl-2H-[1,2,3]triazolo[4,5-b]pyridine'),
    ('Cc1ccnc2n[nH]nc12',   '7-methyl-2H-[1,2,3]triazolo[4,5-b]pyridine'),
    # 1H-[1,2,3]triazolo[4,5-b]pyridine
    ('c1cnc2nn[nH]c2c1',    '1H-[1,2,3]triazolo[4,5-b]pyridine'),
    ('Cc1ccc2[nH]nnc2n1',   '5-methyl-1H-[1,2,3]triazolo[4,5-b]pyridine'),
    ('Cc1cnc2nn[nH]c2c1',   '6-methyl-1H-[1,2,3]triazolo[4,5-b]pyridine'),
    ('Cc1ccnc2nn[nH]c12',   '7-methyl-1H-[1,2,3]triazolo[4,5-b]pyridine'),
])
def test_phase617_nh_bicyclics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
