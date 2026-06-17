"""Phase 569: Substituted benzo[a]pyrene naming.
Benzo[a]pyrene: no symmetry, 20 atoms, 12 unique CH positions (1-12).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34",   "benzo[a]pyrene"),
    ("Cc1ccc2ccc3cc4ccccc4c4ccc1c2c34",    "1-methylbenzo[a]pyrene"),
    ("Cc1cc2ccc3cc4ccccc4c4ccc(c1)c2c34",  "2-methylbenzo[a]pyrene"),
    ("Cc1ccc2ccc3c4ccccc4cc4ccc1c2c43",    "3-methylbenzo[a]pyrene"),
    ("Cc1cc2cc3ccccc3c3ccc4cccc1c4c23",    "4-methylbenzo[a]pyrene"),
    ("Cc1cc2cccc3ccc4c5ccccc5cc1c4c32",    "5-methylbenzo[a]pyrene"),
    ("Cc1c2ccccc2c2ccc3cccc4ccc1c2c43",    "6-methylbenzo[a]pyrene"),
    ("Cc1cccc2c1cc1ccc3cccc4ccc2c1c34",    "7-methylbenzo[a]pyrene"),
    ("Cc1ccc2c(c1)cc1ccc3cccc4ccc2c1c34",  "8-methylbenzo[a]pyrene"),
    ("Cc1ccc2cc3ccc4cccc5ccc(c2c1)c3c45",  "9-methylbenzo[a]pyrene"),
    ("Cc1cccc2cc3ccc4cccc5ccc(c12)c3c45",  "10-methylbenzo[a]pyrene"),
    ("Cc1cc2cccc3ccc4cc5ccccc5c1c4c32",    "11-methylbenzo[a]pyrene"),
    ("Cc1cc2c3ccccc3cc3ccc4cccc1c4c32",    "12-methylbenzo[a]pyrene"),
])
def test_phase569_benzo_a_pyrene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
