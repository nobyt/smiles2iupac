"""Phase 564: Substituted benzo[f]quinoline, benzo[h]quinoline,
benzo[f]isoquinoline, and benzo[h]isoquinoline naming.
All are 14-atom phenanthrene-skeleton 6+6+6 tricyclics with N in one outer ring.
benzo[f]quinoline: N at position 5; substitutable C positions 1-4,6-10.
benzo[h]quinoline: N at position 10; substitutable C positions 1-9.
benzo[f]isoquinoline: N at position 6; substitutable C positions 1-5,7-10.
benzo[h]isoquinoline: N at position 7; substitutable C positions 1-6,8-10.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[f]quinoline (N at pos 5; C positions 1-4,6-10)
    ("c1ccc2c(c1)ccc1ncccc12",          "benzo[f]quinoline"),
    ("Cc1cccc2ccc3ncccc3c12",           "1-methylbenzo[f]quinoline"),
    ("Cc1ccnc2ccc3ccccc3c12",           "2-methylbenzo[f]quinoline"),
    ("Cc1cnc2ccc3ccccc3c2c1",           "3-methylbenzo[f]quinoline"),
    ("Cc1ccc2c(ccc3ccccc32)n1",         "4-methylbenzo[f]quinoline"),
    ("Cc1cc2ccccc2c2cccnc12",           "6-methylbenzo[f]quinoline"),
    ("Cc1cc2ncccc2c2ccccc12",           "7-methylbenzo[f]quinoline"),
    ("Cc1cccc2c1ccc1ncccc12",           "8-methylbenzo[f]quinoline"),
    ("Cc1ccc2c(ccc3ncccc32)c1",         "9-methylbenzo[f]quinoline"),
    ("Cc1ccc2ccc3ncccc3c2c1",           "10-methylbenzo[f]quinoline"),
    # benzo[h]quinoline (N at pos 10; C positions 1-9)
    ("c1ccc2c(c1)ccc1cccnc12",          "benzo[h]quinoline"),
    ("Cc1cccc2ccc3cccnc3c12",           "1-methylbenzo[h]quinoline"),
    ("Cc1ccc2ccc3cccnc3c2c1",           "2-methylbenzo[h]quinoline"),
    ("Cc1ccc2c(ccc3cccnc32)c1",         "3-methylbenzo[h]quinoline"),
    ("Cc1cccc2c1ccc1cccnc12",           "4-methylbenzo[h]quinoline"),
    ("Cc1cc2cccnc2c2ccccc12",           "5-methylbenzo[h]quinoline"),
    ("Cc1cc2ccccc2c2ncccc12",           "6-methylbenzo[h]quinoline"),
    ("Cc1ccnc2c1ccc1ccccc12",           "7-methylbenzo[h]quinoline"),
    ("Cc1cnc2c(ccc3ccccc32)c1",         "8-methylbenzo[h]quinoline"),
    ("Cc1ccc2ccc3ccccc3c2n1",           "9-methylbenzo[h]quinoline"),
    # benzo[f]isoquinoline (N at pos 6; C positions 1-5,7-10)
    ("c1ccc2c(c1)ccc1cnccc12",          "benzo[f]isoquinoline"),
    ("Cc1ccc2c(ccc3cnccc32)c1",         "1-methylbenzo[f]isoquinoline"),
    ("Cc1ccc2ccc3cnccc3c2c1",           "2-methylbenzo[f]isoquinoline"),
    ("Cc1cccc2ccc3cnccc3c12",           "3-methylbenzo[f]isoquinoline"),
    ("Cc1cncc2ccc3ccccc3c12",           "4-methylbenzo[f]isoquinoline"),
    ("Cc1cc2c(ccc3ccccc32)cn1",         "5-methylbenzo[f]isoquinoline"),
    ("Cc1nccc2c1ccc1ccccc12",           "7-methylbenzo[f]isoquinoline"),
    ("Cc1cc2ccccc2c2ccncc12",           "8-methylbenzo[f]isoquinoline"),
    ("Cc1cc2cnccc2c2ccccc12",           "9-methylbenzo[f]isoquinoline"),
    ("Cc1cccc2c1ccc1cnccc12",           "10-methylbenzo[f]isoquinoline"),
    # benzo[h]isoquinoline (N at pos 7; C positions 1-6,8-10)
    ("c1ccc2c(c1)ccc1ccncc12",          "benzo[h]isoquinoline"),
    ("Cc1cc2ccncc2c2ccccc12",           "1-methylbenzo[h]isoquinoline"),
    ("Cc1cccc2c1ccc1ccncc12",           "2-methylbenzo[h]isoquinoline"),
    ("Cc1ccc2c(ccc3ccncc32)c1",         "3-methylbenzo[h]isoquinoline"),
    ("Cc1ccc2ccc3ccncc3c2c1",           "4-methylbenzo[h]isoquinoline"),
    ("Cc1cccc2ccc3ccncc3c12",           "5-methylbenzo[h]isoquinoline"),
    ("Cc1nccc2ccc3ccccc3c12",           "6-methylbenzo[h]isoquinoline"),
    ("Cc1cc2ccc3ccccc3c2cn1",           "8-methylbenzo[h]isoquinoline"),
    ("Cc1cncc2c1ccc1ccccc12",           "9-methylbenzo[h]isoquinoline"),
    ("Cc1cc2ccccc2c2cnccc12",           "10-methylbenzo[h]isoquinoline"),
])
def test_phase564_benzo_quinolines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
