"""Phase 430: benzo[f]isoquinoline and benzo[g]isoquinoline retained names (IUPAC 2013 P-31.1.3).

Two C13H9N tricyclic systems — linear/angular fusion of benzene with isoquinoline.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[f]isoquinoline — angular, N in peripheral ring adjacent to central junction
    ("c1ccc2c(c1)ccc1cnccc12",      "benzo[f]isoquinoline"),
    # benzo[g]isoquinoline — linear, N in terminal ring
    ("c1ccc2cc3cnccc3cc2c1",        "benzo[g]isoquinoline"),
    # regression: benzo[f]quinoline unchanged (Phase 429)
    ("c1ccc2c(c1)ccc1ncccc12",      "benzo[f]quinoline"),
    # regression: benzo[h]quinoline unchanged (Phase 429)
    ("c1ccc2c(c1)ccc1cccnc12",      "benzo[h]quinoline"),
    # regression: isoquinoline unchanged
    ("c1ccc2cnccc2c1",              "isoquinoline"),
    # regression: acridine unchanged
    ("c1ccc2nc3ccccc3cc2c1",        "acridine"),
    # regression: naphthalene unchanged
    ("c1ccc2ccccc2c1",              "naphthalene"),
])
def test_phase430_benzoisoquinoline(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
