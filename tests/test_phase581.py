"""Phase 581: Substituted 1H-naphtho[1,2-d][1,2,3]triazole, naphtho[2,3-d][1,2,3]oxadiazole,
naphtho[1,2-d][1,2,3]oxadiazole, naphtho[2,1-d][1,2,3]oxadiazole,
naphtho[2,3-d][1,2,3]thiadiazole, naphtho[1,2-d][1,2,3]thiadiazole,
and naphtho[2,1-d][1,2,3]thiadiazole naming.
Heteroatoms: N/O/S (all at None); substitutable C positions: 4-9 (6 positions).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-naphtho[1,2-d][1,2,3]triazole
    ("c1ccc2c(c1)ccc1nn[nH]c12",        "1H-naphtho[1,2-d][1,2,3]triazole"),
    ("Cc1cc2ccccc2c2[nH]nnc12",         "4-methyl-1H-naphtho[1,2-d][1,2,3]triazole"),
    ("Cc1cc2nn[nH]c2c2ccccc12",         "5-methyl-1H-naphtho[1,2-d][1,2,3]triazole"),
    ("Cc1cccc2c1ccc1nn[nH]c12",         "6-methyl-1H-naphtho[1,2-d][1,2,3]triazole"),
    ("Cc1ccc2c(ccc3nn[nH]c32)c1",       "7-methyl-1H-naphtho[1,2-d][1,2,3]triazole"),
    ("Cc1ccc2ccc3nn[nH]c3c2c1",         "8-methyl-1H-naphtho[1,2-d][1,2,3]triazole"),
    ("Cc1cccc2ccc3nn[nH]c3c12",         "9-methyl-1H-naphtho[1,2-d][1,2,3]triazole"),
    # naphtho[2,3-d][1,2,3]oxadiazole
    ("c1ccc2cc3onnc3cc2c1",             "naphtho[2,3-d][1,2,3]oxadiazole"),
    ("Cc1c2ccccc2cc2onnc12",            "4-methylnaphtho[2,3-d][1,2,3]oxadiazole"),
    ("Cc1cccc2cc3onnc3cc12",            "5-methylnaphtho[2,3-d][1,2,3]oxadiazole"),
    ("Cc1ccc2cc3onnc3cc2c1",            "6-methylnaphtho[2,3-d][1,2,3]oxadiazole"),
    ("Cc1ccc2cc3nnoc3cc2c1",            "7-methylnaphtho[2,3-d][1,2,3]oxadiazole"),
    ("Cc1cccc2cc3nnoc3cc12",            "8-methylnaphtho[2,3-d][1,2,3]oxadiazole"),
    ("Cc1c2ccccc2cc2nnoc12",            "9-methylnaphtho[2,3-d][1,2,3]oxadiazole"),
    # naphtho[1,2-d][1,2,3]oxadiazole
    ("c1ccc2c(c1)ccc1onnc12",           "naphtho[1,2-d][1,2,3]oxadiazole"),
    ("Cc1cc2ccccc2c2nnoc12",            "4-methylnaphtho[1,2-d][1,2,3]oxadiazole"),
    ("Cc1cc2onnc2c2ccccc12",            "5-methylnaphtho[1,2-d][1,2,3]oxadiazole"),
    ("Cc1cccc2c1ccc1onnc12",            "6-methylnaphtho[1,2-d][1,2,3]oxadiazole"),
    ("Cc1ccc2c(ccc3onnc32)c1",          "7-methylnaphtho[1,2-d][1,2,3]oxadiazole"),
    ("Cc1ccc2ccc3onnc3c2c1",            "8-methylnaphtho[1,2-d][1,2,3]oxadiazole"),
    ("Cc1cccc2ccc3onnc3c12",            "9-methylnaphtho[1,2-d][1,2,3]oxadiazole"),
    # naphtho[2,1-d][1,2,3]oxadiazole
    ("c1ccc2c(c1)ccc1nnoc12",           "naphtho[2,1-d][1,2,3]oxadiazole"),
    ("Cc1cc2ccccc2c2onnc12",            "4-methylnaphtho[2,1-d][1,2,3]oxadiazole"),
    ("Cc1cc2nnoc2c2ccccc12",            "5-methylnaphtho[2,1-d][1,2,3]oxadiazole"),
    ("Cc1cccc2c1ccc1nnoc12",            "6-methylnaphtho[2,1-d][1,2,3]oxadiazole"),
    ("Cc1ccc2c(ccc3nnoc32)c1",          "7-methylnaphtho[2,1-d][1,2,3]oxadiazole"),
    ("Cc1ccc2ccc3nnoc3c2c1",            "8-methylnaphtho[2,1-d][1,2,3]oxadiazole"),
    ("Cc1cccc2ccc3nnoc3c12",            "9-methylnaphtho[2,1-d][1,2,3]oxadiazole"),
    # naphtho[2,3-d][1,2,3]thiadiazole
    ("c1ccc2cc3snnc3cc2c1",             "naphtho[2,3-d][1,2,3]thiadiazole"),
    ("Cc1c2ccccc2cc2snnc12",            "4-methylnaphtho[2,3-d][1,2,3]thiadiazole"),
    ("Cc1cccc2cc3snnc3cc12",            "5-methylnaphtho[2,3-d][1,2,3]thiadiazole"),
    ("Cc1ccc2cc3snnc3cc2c1",            "6-methylnaphtho[2,3-d][1,2,3]thiadiazole"),
    ("Cc1ccc2cc3nnsc3cc2c1",            "7-methylnaphtho[2,3-d][1,2,3]thiadiazole"),
    ("Cc1cccc2cc3nnsc3cc12",            "8-methylnaphtho[2,3-d][1,2,3]thiadiazole"),
    ("Cc1c2ccccc2cc2nnsc12",            "9-methylnaphtho[2,3-d][1,2,3]thiadiazole"),
    # naphtho[1,2-d][1,2,3]thiadiazole
    ("c1ccc2c(c1)ccc1snnc12",           "naphtho[1,2-d][1,2,3]thiadiazole"),
    ("Cc1cc2ccccc2c2nnsc12",            "4-methylnaphtho[1,2-d][1,2,3]thiadiazole"),
    ("Cc1cc2snnc2c2ccccc12",            "5-methylnaphtho[1,2-d][1,2,3]thiadiazole"),
    ("Cc1cccc2c1ccc1snnc12",            "6-methylnaphtho[1,2-d][1,2,3]thiadiazole"),
    ("Cc1ccc2c(ccc3snnc32)c1",          "7-methylnaphtho[1,2-d][1,2,3]thiadiazole"),
    ("Cc1ccc2ccc3snnc3c2c1",            "8-methylnaphtho[1,2-d][1,2,3]thiadiazole"),
    ("Cc1cccc2ccc3snnc3c12",            "9-methylnaphtho[1,2-d][1,2,3]thiadiazole"),
    # naphtho[2,1-d][1,2,3]thiadiazole
    ("c1ccc2c(c1)ccc1nnsc12",           "naphtho[2,1-d][1,2,3]thiadiazole"),
    ("Cc1cc2ccccc2c2snnc12",            "4-methylnaphtho[2,1-d][1,2,3]thiadiazole"),
    ("Cc1cc2nnsc2c2ccccc12",            "5-methylnaphtho[2,1-d][1,2,3]thiadiazole"),
    ("Cc1cccc2c1ccc1nnsc12",            "6-methylnaphtho[2,1-d][1,2,3]thiadiazole"),
    ("Cc1ccc2c(ccc3nnsc32)c1",          "7-methylnaphtho[2,1-d][1,2,3]thiadiazole"),
    ("Cc1ccc2ccc3nnsc3c2c1",            "8-methylnaphtho[2,1-d][1,2,3]thiadiazole"),
    ("Cc1cccc2ccc3nnsc3c12",            "9-methylnaphtho[2,1-d][1,2,3]thiadiazole"),
])
def test_phase581_naphtho_triazole_oxadiazole_thiadiazole(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
