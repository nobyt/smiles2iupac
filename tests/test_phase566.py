"""Phase 566: Substituted benzo[c]cinnoline, benzo[f]cinnoline, benzo[h]cinnoline,
benzo[f]phthalazine, benzo[f]quinoxaline, benzo[g]cinnoline, benzo[g]phthalazine,
and benzo[g]quinoxaline naming.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[c]cinnoline (C2-symmetric; 4 unique C positions: 1,2,3,4)
    ("c1ccc2c(c1)nnc1ccccc12",          "benzo[c]cinnoline"),
    ("Cc1cccc2nnc3ccccc3c12",           "1-methylbenzo[c]cinnoline"),
    ("Cc1ccc2nnc3ccccc3c2c1",           "2-methylbenzo[c]cinnoline"),
    ("Cc1ccc2c(c1)nnc1ccccc12",         "3-methylbenzo[c]cinnoline"),
    ("Cc1cccc2c1nnc1ccccc12",           "4-methylbenzo[c]cinnoline"),
    # benzo[f]cinnoline (8 unique C positions: 1,2,3,6,7,8,9,10)
    ("c1ccc2c(c1)ccc1nnccc12",          "benzo[f]cinnoline"),
    ("Cc1cccc2ccc3nnccc3c12",           "1-methylbenzo[f]cinnoline"),
    ("Cc1cnnc2ccc3ccccc3c12",           "2-methylbenzo[f]cinnoline"),
    ("Cc1cc2c(ccc3ccccc32)nn1",         "3-methylbenzo[f]cinnoline"),
    ("Cc1cc2ccccc2c2ccnnc12",           "6-methylbenzo[f]cinnoline"),
    ("Cc1cc2nnccc2c2ccccc12",           "7-methylbenzo[f]cinnoline"),
    ("Cc1cccc2c1ccc1nnccc12",           "8-methylbenzo[f]cinnoline"),
    ("Cc1ccc2c(ccc3nnccc32)c1",         "9-methylbenzo[f]cinnoline"),
    ("Cc1ccc2ccc3nnccc3c2c1",           "10-methylbenzo[f]cinnoline"),
    # benzo[h]cinnoline (8 unique C positions: 1,2,3,4,5,6,7,8)
    ("c1ccc2c(c1)ccc1ccnnc12",          "benzo[h]cinnoline"),
    ("Cc1cccc2ccc3ccnnc3c12",           "1-methylbenzo[h]cinnoline"),
    ("Cc1ccc2ccc3ccnnc3c2c1",           "2-methylbenzo[h]cinnoline"),
    ("Cc1ccc2c(ccc3ccnnc32)c1",         "3-methylbenzo[h]cinnoline"),
    ("Cc1cccc2c1ccc1ccnnc12",           "4-methylbenzo[h]cinnoline"),
    ("Cc1cc2ccnnc2c2ccccc12",           "5-methylbenzo[h]cinnoline"),
    ("Cc1cc2ccccc2c2nnccc12",           "6-methylbenzo[h]cinnoline"),
    ("Cc1cnnc2c1ccc1ccccc12",           "7-methylbenzo[h]cinnoline"),
    ("Cc1cc2ccc3ccccc3c2nn1",           "8-methylbenzo[h]cinnoline"),
    # benzo[f]phthalazine (8 unique C positions: 1,2,5,6,7,8,9,10)
    ("c1ccc2c(c1)ccc1cnncc12",          "benzo[f]phthalazine"),
    ("Cc1cccc2ccc3cnncc3c12",           "1-methylbenzo[f]phthalazine"),
    ("Cc1nncc2ccc3ccccc3c12",           "2-methylbenzo[f]phthalazine"),
    ("Cc1nncc2c1ccc1ccccc12",           "5-methylbenzo[f]phthalazine"),
    ("Cc1cc2ccccc2c2cnncc12",           "6-methylbenzo[f]phthalazine"),
    ("Cc1cc2cnncc2c2ccccc12",           "7-methylbenzo[f]phthalazine"),
    ("Cc1cccc2c1ccc1cnncc12",           "8-methylbenzo[f]phthalazine"),
    ("Cc1ccc2c(ccc3cnncc32)c1",         "9-methylbenzo[f]phthalazine"),
    ("Cc1ccc2ccc3cnncc3c2c1",           "10-methylbenzo[f]phthalazine"),
    # benzo[f]quinoxaline (8 unique C positions: 1,3,4,6,7,8,9,10)
    ("c1ccc2c(c1)ccc1nccnc12",          "benzo[f]quinoxaline"),
    ("Cc1cccc2ccc3nccnc3c12",           "1-methylbenzo[f]quinoxaline"),
    ("Cc1cnc2ccc3ccccc3c2n1",           "3-methylbenzo[f]quinoxaline"),
    ("Cc1cnc2c(ccc3ccccc32)n1",         "4-methylbenzo[f]quinoxaline"),
    ("Cc1cc2ccccc2c2nccnc12",           "6-methylbenzo[f]quinoxaline"),
    ("Cc1cc2nccnc2c2ccccc12",           "7-methylbenzo[f]quinoxaline"),
    ("Cc1cccc2c1ccc1nccnc12",           "8-methylbenzo[f]quinoxaline"),
    ("Cc1ccc2c(ccc3nccnc32)c1",         "9-methylbenzo[f]quinoxaline"),
    ("Cc1ccc2ccc3nccnc3c2c1",           "10-methylbenzo[f]quinoxaline"),
    # benzo[g]cinnoline (8 unique C positions: 1,2,3,4,5,6,9,10)
    ("c1ccc2cc3nnccc3cc2c1",            "benzo[g]cinnoline"),
    ("Cc1cccc2cc3ccnnc3cc12",           "1-methylbenzo[g]cinnoline"),
    ("Cc1ccc2cc3ccnnc3cc2c1",           "2-methylbenzo[g]cinnoline"),
    ("Cc1ccc2cc3nnccc3cc2c1",           "3-methylbenzo[g]cinnoline"),
    ("Cc1cccc2cc3nnccc3cc12",           "4-methylbenzo[g]cinnoline"),
    ("Cc1cnnc2cc3ccccc3cc12",           "5-methylbenzo[g]cinnoline"),
    ("Cc1cc2cc3ccccc3cc2nn1",           "6-methylbenzo[g]cinnoline"),
    ("Cc1c2ccccc2cc2nnccc12",           "9-methylbenzo[g]cinnoline"),
    ("Cc1c2ccccc2cc2ccnnc12",           "10-methylbenzo[g]cinnoline"),
    # benzo[g]phthalazine (C2-symmetric; 4 unique C positions: 1,2,5,9)
    ("c1ccc2cc3cnncc3cc2c1",            "benzo[g]phthalazine"),
    ("Cc1cccc2cc3cnncc3cc12",           "1-methylbenzo[g]phthalazine"),
    ("Cc1ccc2cc3cnncc3cc2c1",           "2-methylbenzo[g]phthalazine"),
    ("Cc1nncc2cc3ccccc3cc12",           "5-methylbenzo[g]phthalazine"),
    ("Cc1c2ccccc2cc2cnncc12",           "9-methylbenzo[g]phthalazine"),
    # benzo[g]quinoxaline (C2-symmetric; 4 unique C positions: 1,2,6,9)
    ("c1ccc2cc3nccnc3cc2c1",            "benzo[g]quinoxaline"),
    ("Cc1cccc2cc3nccnc3cc12",           "1-methylbenzo[g]quinoxaline"),
    ("Cc1ccc2cc3nccnc3cc2c1",           "2-methylbenzo[g]quinoxaline"),
    ("Cc1cnc2cc3ccccc3cc2n1",           "6-methylbenzo[g]quinoxaline"),
    ("Cc1c2ccccc2cc2nccnc12",           "9-methylbenzo[g]quinoxaline"),
])
def test_phase566_benzo_diazines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
