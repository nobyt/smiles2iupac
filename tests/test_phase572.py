"""Phase 572: Substituted benzo[g]isoquinoline naming.
Benzo[g]isoquinoline: N at position 2, 14 atoms, 9 unique C positions (1,3-10).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1ccc2cc3cnccc3cc2c1",           "benzo[g]isoquinoline"),
    ("Cc1nccc2cc3ccccc3cc12",          "1-methylbenzo[g]isoquinoline"),
    ("Cc1cc2cc3ccccc3cc2cn1",          "3-methylbenzo[g]isoquinoline"),
    ("Cc1cncc2cc3ccccc3cc12",          "4-methylbenzo[g]isoquinoline"),
    ("Cc1c2ccccc2cc2cnccc12",          "5-methylbenzo[g]isoquinoline"),
    ("Cc1cccc2cc3cnccc3cc12",          "6-methylbenzo[g]isoquinoline"),
    ("Cc1ccc2cc3cnccc3cc2c1",          "7-methylbenzo[g]isoquinoline"),
    ("Cc1ccc2cc3ccncc3cc2c1",          "8-methylbenzo[g]isoquinoline"),
    ("Cc1cccc2cc3ccncc3cc12",          "9-methylbenzo[g]isoquinoline"),
    ("Cc1c2ccccc2cc2ccncc12",          "10-methylbenzo[g]isoquinoline"),
])
def test_phase572_benzo_g_isoquinoline(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
