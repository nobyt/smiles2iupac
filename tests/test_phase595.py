"""Phase 595: Substituted indolizine, imidazo[1,2-a]pyridine, imidazo[1,5-a]pyridine,
pyrazolo[1,5-a]pyridine, [1,2,4]triazolo[4,3-a]pyridine,
[1,2,4]triazolo[1,5-a]pyridine, [1,2,3]triazolo[1,5-a]pyridine,
and tetrazolo[1,5-a]pyridine naming.
N-bridged 5+6 bicyclics.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # indolizine (N at 4; sub C: 1,2,3,5,6,7,8)
    ("c1ccn2cccc2c1",          "indolizine"),
    ("Cc1ccn2ccccc12",         "1-methylindolizine"),
    ("Cc1cc2ccccn2c1",         "2-methylindolizine"),
    ("Cc1ccc2ccccn12",         "3-methylindolizine"),
    ("Cc1cccc2cccn12",         "5-methylindolizine"),
    ("Cc1ccc2cccn2c1",         "6-methylindolizine"),
    ("Cc1ccn2cccc2c1",         "7-methylindolizine"),
    ("Cc1cccn2cccc12",         "8-methylindolizine"),
    # imidazo[1,2-a]pyridine (N at 4,9; sub C: 2,3,5,6,7,8)
    ("c1ccn2ccnc2c1",          "imidazo[1,2-a]pyridine"),
    ("Cc1cn2ccccc2n1",         "2-methylimidazo[1,2-a]pyridine"),
    ("Cc1cnc2ccccn12",         "3-methylimidazo[1,2-a]pyridine"),
    ("Cc1cccc2nccn12",         "5-methylimidazo[1,2-a]pyridine"),
    ("Cc1ccc2nccn2c1",         "6-methylimidazo[1,2-a]pyridine"),
    ("Cc1ccn2ccnc2c1",         "7-methylimidazo[1,2-a]pyridine"),
    ("Cc1cccn2ccnc12",         "8-methylimidazo[1,2-a]pyridine"),
    # imidazo[1,5-a]pyridine (N at 4,5,9; sub C: 1,3,5,6,7,8)
    ("c1ccn2cncc2c1",          "imidazo[1,5-a]pyridine"),
    ("Cc1ncn2ccccc12",         "1-methylimidazo[1,5-a]pyridine"),
    ("Cc1ncc2ccccn12",         "3-methylimidazo[1,5-a]pyridine"),
    ("Cc1cccc2cncn12",         "5-methylimidazo[1,5-a]pyridine"),
    ("Cc1ccc2cncn2c1",         "6-methylimidazo[1,5-a]pyridine"),
    ("Cc1ccn2cncc2c1",         "7-methylimidazo[1,5-a]pyridine"),
    ("Cc1cccn2cncc12",         "8-methylimidazo[1,5-a]pyridine"),
    # pyrazolo[1,5-a]pyridine (N at 4,5; sub C: 2,3,5,6,7)
    ("c1ccn2nccc2c1",          "pyrazolo[1,5-a]pyridine"),
    ("Cc1cc2ccccn2n1",         "2-methylpyrazolo[1,5-a]pyridine"),
    ("Cc1cnn2ccccc12",         "3-methylpyrazolo[1,5-a]pyridine"),
    ("Cc1ccn2nccc2c1",         "5-methylpyrazolo[1,5-a]pyridine"),
    ("Cc1ccc2ccnn2c1",         "6-methylpyrazolo[1,5-a]pyridine"),
    ("Cc1cccc2ccnn12",         "7-methylpyrazolo[1,5-a]pyridine"),
    # [1,2,4]triazolo[4,3-a]pyridine (N at 4,6,7; sub C: 3,5,6,7,8)
    ("c1ccn2cnnc2c1",          "[1,2,4]triazolo[4,3-a]pyridine"),
    ("Cc1nnc2ccccn12",         "3-methyl[1,2,4]triazolo[4,3-a]pyridine"),
    ("Cc1cccc2nncn12",         "5-methyl[1,2,4]triazolo[4,3-a]pyridine"),
    ("Cc1ccc2nncn2c1",         "6-methyl[1,2,4]triazolo[4,3-a]pyridine"),
    ("Cc1ccn2cnnc2c1",         "7-methyl[1,2,4]triazolo[4,3-a]pyridine"),
    ("Cc1cccn2cnnc12",         "8-methyl[1,2,4]triazolo[4,3-a]pyridine"),
    # [1,2,4]triazolo[1,5-a]pyridine (N at 4,5,7; sub C: 2,5,6,7,8)
    ("c1ccn2ncnc2c1",          "[1,2,4]triazolo[1,5-a]pyridine"),
    ("Cc1nc2ccccn2n1",         "2-methyl[1,2,4]triazolo[1,5-a]pyridine"),
    ("Cc1cccc2ncnn12",         "5-methyl[1,2,4]triazolo[1,5-a]pyridine"),
    ("Cc1ccc2ncnn2c1",         "6-methyl[1,2,4]triazolo[1,5-a]pyridine"),
    ("Cc1ccn2ncnc2c1",         "7-methyl[1,2,4]triazolo[1,5-a]pyridine"),
    ("Cc1cccn2ncnc12",         "8-methyl[1,2,4]triazolo[1,5-a]pyridine"),
    # [1,2,3]triazolo[1,5-a]pyridine (N at 4,5,6; sub C: 3,5,6,7)
    ("c1ccn2nncc2c1",          "[1,2,3]triazolo[1,5-a]pyridine"),
    ("Cc1nnn2ccccc12",         "3-methyl[1,2,3]triazolo[1,5-a]pyridine"),
    ("Cc1ccn2nncc2c1",         "5-methyl[1,2,3]triazolo[1,5-a]pyridine"),
    ("Cc1ccc2cnnn2c1",         "6-methyl[1,2,3]triazolo[1,5-a]pyridine"),
    ("Cc1cccc2cnnn12",         "7-methyl[1,2,3]triazolo[1,5-a]pyridine"),
    # tetrazolo[1,5-a]pyridine (N at 4,5,6,7; sub C: 5,6,7,8)
    ("c1ccn2nnnc2c1",          "tetrazolo[1,5-a]pyridine"),
    ("Cc1cccc2nnnn12",         "5-methyltetrazolo[1,5-a]pyridine"),
    ("Cc1ccc2nnnn2c1",         "6-methyltetrazolo[1,5-a]pyridine"),
    ("Cc1ccn2nnnc2c1",         "7-methyltetrazolo[1,5-a]pyridine"),
    ("Cc1cccn2nnnc12",         "8-methyltetrazolo[1,5-a]pyridine"),
])
def test_phase595_n_bridged_bicyclics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
