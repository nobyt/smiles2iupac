"""Phase 613: monocyclic heteroaromatics — triazines, tetrazine, isothiazole, isoxazole, oxazole, thiazole, thieno isomers"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2,3-triazine
    ('c1cnnnc1', '1,2,3-triazine'),
    ('Cc1ccnnn1', '4-methyl-1,2,3-triazine'),
    ('Cc1cnnnc1', '5-methyl-1,2,3-triazine'),
    # 1,2,4,5-tetrazine
    ('c1nncnn1', '1,2,4,5-tetrazine'),
    ('Cc1nncnn1', '3-methyl-1,2,4,5-tetrazine'),
    # 1,2,4-triazine
    ('c1cnncn1', '1,2,4-triazine'),
    ('Cc1nccnn1', '3-methyl-1,2,4-triazine'),
    ('Cc1cnncn1', '5-methyl-1,2,4-triazine'),
    ('Cc1cncnn1', '6-methyl-1,2,4-triazine'),
    # 1,3,5-triazine
    ('c1ncncn1', '1,3,5-triazine'),
    ('Cc1ncncn1', '2-methyl-1,3,5-triazine'),
    # isothiazole
    ('c1cnsc1', 'isothiazole'),
    ('Cc1ccsn1', '3-methylisothiazole'),
    ('Cc1cnsc1', '4-methylisothiazole'),
    ('Cc1ccns1', '5-methylisothiazole'),
    # isoxazole
    ('c1cnoc1', 'isoxazole'),
    ('Cc1ccon1', '3-methylisoxazole'),
    ('Cc1cnoc1', '4-methylisoxazole'),
    ('Cc1ccno1', '5-methylisoxazole'),
    # oxazole
    ('c1cocn1', 'oxazole'),
    ('Cc1ncco1', '2-methyloxazole'),
    ('Cc1cocn1', '4-methyloxazole'),
    ('Cc1cnco1', '5-methyloxazole'),
    # thiazole
    ('c1cscn1', 'thiazole'),
    ('Cc1nccs1', '2-methylthiazole'),
    ('Cc1cscn1', '4-methylthiazole'),
    ('Cc1cncs1', '5-methylthiazole'),
    # thieno[2,3-b]thiophene
    ('c1cc2ccsc2s1', 'thieno[2,3-b]thiophene'),
    ('Cc1cc2ccsc2s1', '2-methylthieno[2,3-b]thiophene'),
    ('Cc1csc2sccc12', '3-methylthieno[2,3-b]thiophene'),
    # thieno[3,2-b]thiophene
    ('c1cc2sccc2s1', 'thieno[3,2-b]thiophene'),
    ('Cc1cc2sccc2s1', '2-methylthieno[3,2-b]thiophene'),
    ('Cc1csc2ccsc12', '3-methylthieno[3,2-b]thiophene'),
    # thieno[3,2-b]pyridine
    ('c1cnc2ccsc2c1', 'thieno[3,2-b]pyridine'),
    ('Cc1cc2ncccc2s1', '2-methylthieno[3,2-b]pyridine'),
    ('Cc1csc2cccnc12', '3-methylthieno[3,2-b]pyridine'),
    ('Cc1ccc2sccc2n1', '5-methylthieno[3,2-b]pyridine'),
    ('Cc1cnc2ccsc2c1', '6-methylthieno[3,2-b]pyridine'),
    ('Cc1ccnc2ccsc12', '7-methylthieno[3,2-b]pyridine'),
    # thieno[2,3-b]pyridine
    ('c1cnc2sccc2c1', 'thieno[2,3-b]pyridine'),
    ('Cc1cc2cccnc2s1', '2-methylthieno[2,3-b]pyridine'),
    ('Cc1csc2ncccc12', '3-methylthieno[2,3-b]pyridine'),
    ('Cc1ccnc2sccc12', '4-methylthieno[2,3-b]pyridine'),
    ('Cc1cnc2sccc2c1', '5-methylthieno[2,3-b]pyridine'),
    ('Cc1ccc2ccsc2n1', '6-methylthieno[2,3-b]pyridine'),
])
def test_phase613_heteroaromatics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
