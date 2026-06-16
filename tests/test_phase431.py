"""Phase 431: benz[a]acridine and benz[c]acridine retained names (IUPAC 2013 P-31.1.3).

Two C17H11N tetracyclic systems formed by fusion of benzene with acridine.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benz[a]acridine — benzo fused at [a] bond of acridine (N-adjacent ring)
    ("c1ccc2nc3ccc4ccccc4c3cc2c1",    "benz[a]acridine"),
    # benz[c]acridine — benzo fused at [c] bond of acridine
    ("c1ccc2cc3nc4ccccc4cc3cc2c1",    "benz[c]acridine"),
    # regression: acridine unchanged
    ("c1ccc2nc3ccccc3cc2c1",           "acridine"),
    # regression: benzo[f]quinoline unchanged (Phase 429)
    ("c1ccc2c(c1)ccc1ncccc12",         "benzo[f]quinoline"),
    # regression: anthracene unchanged (tetracyclic without N)
    ("c1ccc2cc3ccccc3cc2c1",           "anthracene"),
])
def test_phase431_benzacridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
