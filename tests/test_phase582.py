"""Phase 582: Substituted benzo[g]quinoline naming.
benzo[g]quinoline = naphtho[2,3-b]pyridine (IUPAC preferred retained name).
Substitutable C positions: 2-10 (9 positions).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[g]quinoline
    ("c1ccc2cc3ncccc3cc2c1",        "benzo[g]quinoline"),
    ("Cc1ccc2cc3ccccc3cc2n1",       "2-methylbenzo[g]quinoline"),
    ("Cc1cnc2cc3ccccc3cc2c1",       "3-methylbenzo[g]quinoline"),
    ("Cc1ccnc2cc3ccccc3cc12",       "4-methylbenzo[g]quinoline"),
    ("Cc1c2ccccc2cc2ncccc12",       "5-methylbenzo[g]quinoline"),
    ("Cc1cccc2cc3ncccc3cc12",       "6-methylbenzo[g]quinoline"),
    ("Cc1ccc2cc3ncccc3cc2c1",       "7-methylbenzo[g]quinoline"),
    ("Cc1ccc2cc3cccnc3cc2c1",       "8-methylbenzo[g]quinoline"),
    ("Cc1cccc2cc3cccnc3cc12",       "9-methylbenzo[g]quinoline"),
    ("Cc1c2ccccc2cc2cccnc12",       "10-methylbenzo[g]quinoline"),
])
def test_phase582_benzo_g_quinoline(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
