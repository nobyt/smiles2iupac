"""Phase 429: benzo[f]quinoline and benzo[h]quinoline retained names (IUPAC 2013 P-31.1.3).

Two C13H9N tricyclic systems formed by linear fusion of benzene with quinoline.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[f]quinoline — angular tricycle, N adjacent to the central junction
    ("c1ccc2c(c1)ccc1ncccc12",       "benzo[f]quinoline"),
    # benzo[h]quinoline — angular tricycle, N remote from the central junction
    ("c1ccc2c(c1)ccc1cccnc12",       "benzo[h]quinoline"),
    # regression: acridine unchanged (N in central ring)
    ("c1ccc2nc3ccccc3cc2c1",          "acridine"),
    # regression: phenanthridine unchanged (Phase 134)
    ("c1ccc2c(c1)cnc1ccccc12",        "phenanthridine"),
    # regression: quinoline unchanged
    ("c1ccc2ncccc2c1",                "quinoline"),
    # regression: naphthalene unchanged
    ("c1ccc2ccccc2c1",                "naphthalene"),
])
def test_phase429_benzoquinoline(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
