"""Phase 567: Substituted pyrene, triphenylene, and tetracene naming.
Pyrene: D2h-symmetric, 3 unique C environments: {1,3,6,8}, {2,7}, {4,5,9,10}.
Triphenylene: D3h-symmetric, 2 unique C environments: {1,5,9}, {2,6,10}.
Tetracene: D2h-symmetric, 3 unique C environments: {1,4,7,10}, {2,3,8,9}, {5,6,11,12}.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrene (D2h; positions 1,2,3,4,5,6,7,8,9,10)
    ("c1cc2ccc3cccc4ccc(c1)c2c34",       "pyrene"),
    ("Cc1ccc2ccc3cccc4ccc1c2c34",        "1-methylpyrene"),
    ("Cc1cc2ccc3cccc4ccc(c1)c2c34",      "2-methylpyrene"),
    ("Cc1cc2cccc3ccc4cccc1c4c32",        "4-methylpyrene"),
    # triphenylene (D3h; positions 1,2,3,4,5,6,7,8,9,10,11,12)
    ("c1ccc2c(c1)c1ccccc1c1ccccc21",     "triphenylene"),
    ("Cc1cccc2c3ccccc3c3ccccc3c12",      "1-methyltriphenylene"),
    ("Cc1ccc2c3ccccc3c3ccccc3c2c1",      "2-methyltriphenylene"),
    # tetracene (D2h; positions 1,2,3,4,5,6,7,8,9,10,11,12)
    ("c1ccc2cc3cc4ccccc4cc3cc2c1",       "tetracene"),
    ("Cc1cccc2cc3cc4ccccc4cc3cc12",      "1-methyltetracene"),
    ("Cc1ccc2cc3cc4ccccc4cc3cc2c1",      "2-methyltetracene"),
    ("Cc1c2ccccc2cc2cc3ccccc3cc12",      "5-methyltetracene"),
])
def test_phase567_polycyclics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
