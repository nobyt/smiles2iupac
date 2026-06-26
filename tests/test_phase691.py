"""Phase 691: methyl derivatives of benzo[f/g/h]quinoline and benzo[f/g/h]isoquinoline."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[f]quinoline (N at 5; positions 1-4, 6-10 are CH)
    ("c1ccc2c(c1)ccc1ncccc12",        "benzo[f]quinoline"),
    ("Cc1cccc2ccc3ncccc3c12",          "1-methylbenzo[f]quinoline"),
    ("Cc1ccnc2ccc3ccccc3c12",          "2-methylbenzo[f]quinoline"),
    ("Cc1cnc2ccc3ccccc3c2c1",          "3-methylbenzo[f]quinoline"),
    ("Cc1ccc2c(ccc3ccccc32)n1",        "4-methylbenzo[f]quinoline"),
    ("Cc1cc2ccccc2c2cccnc12",          "6-methylbenzo[f]quinoline"),
    ("Cc1cc2ncccc2c2ccccc12",          "7-methylbenzo[f]quinoline"),
    ("Cc1cccc2c1ccc1ncccc12",          "8-methylbenzo[f]quinoline"),
    ("Cc1ccc2c(ccc3ncccc32)c1",        "9-methylbenzo[f]quinoline"),
    ("Cc1ccc2ccc3ncccc3c2c1",          "10-methylbenzo[f]quinoline"),
    # benzo[g]quinoline (N at 1; positions 2-10 are CH)
    ("c1ccc2cc3ncccc3cc2c1",           "benzo[g]quinoline"),
    ("Cc1ccc2cc3ccccc3cc2n1",          "2-methylbenzo[g]quinoline"),
    ("Cc1cnc2cc3ccccc3cc2c1",          "3-methylbenzo[g]quinoline"),
    ("Cc1ccnc2cc3ccccc3cc12",          "4-methylbenzo[g]quinoline"),
    ("Cc1c2ccccc2cc2ncccc12",          "5-methylbenzo[g]quinoline"),
    ("Cc1cccc2cc3ncccc3cc12",          "6-methylbenzo[g]quinoline"),
    ("Cc1ccc2cc3ncccc3cc2c1",          "7-methylbenzo[g]quinoline"),
    ("Cc1ccc2cc3cccnc3cc2c1",          "8-methylbenzo[g]quinoline"),
    ("Cc1cccc2cc3cccnc3cc12",          "9-methylbenzo[g]quinoline"),
    ("Cc1c2ccccc2cc2cccnc12",          "10-methylbenzo[g]quinoline"),
    # benzo[h]quinoline (N at 10; positions 1-9 are CH)
    ("c1ccc2c(c1)ccc1cccnc12",         "benzo[h]quinoline"),
    ("Cc1cccc2ccc3cccnc3c12",          "1-methylbenzo[h]quinoline"),
    ("Cc1ccc2ccc3cccnc3c2c1",          "2-methylbenzo[h]quinoline"),
    ("Cc1ccc2c(ccc3cccnc32)c1",        "3-methylbenzo[h]quinoline"),
    ("Cc1cccc2c1ccc1cccnc12",          "4-methylbenzo[h]quinoline"),
    ("Cc1cc2cccnc2c2ccccc12",          "5-methylbenzo[h]quinoline"),
    ("Cc1cc2ccccc2c2ncccc12",          "6-methylbenzo[h]quinoline"),
    ("Cc1ccnc2c1ccc1ccccc12",          "7-methylbenzo[h]quinoline"),
    ("Cc1cnc2c(ccc3ccccc32)c1",        "8-methylbenzo[h]quinoline"),
    ("Cc1ccc2ccc3ccccc3c2n1",          "9-methylbenzo[h]quinoline"),
    # benzo[f]isoquinoline (N at 6; positions 1-5, 7-10 are CH)
    ("c1ccc2c(c1)ccc1cnccc12",         "benzo[f]isoquinoline"),
    ("Cc1ccc2c(ccc3cnccc32)c1",        "1-methylbenzo[f]isoquinoline"),
    ("Cc1ccc2ccc3cnccc3c2c1",          "2-methylbenzo[f]isoquinoline"),
    ("Cc1cccc2ccc3cnccc3c12",          "3-methylbenzo[f]isoquinoline"),
    ("Cc1cncc2ccc3ccccc3c12",          "4-methylbenzo[f]isoquinoline"),
    ("Cc1cc2c(ccc3ccccc32)cn1",        "5-methylbenzo[f]isoquinoline"),
    ("Cc1nccc2c1ccc1ccccc12",          "7-methylbenzo[f]isoquinoline"),
    ("Cc1cc2ccccc2c2ccncc12",          "8-methylbenzo[f]isoquinoline"),
    ("Cc1cc2cnccc2c2ccccc12",          "9-methylbenzo[f]isoquinoline"),
    ("Cc1cccc2c1ccc1cnccc12",          "10-methylbenzo[f]isoquinoline"),
    # benzo[g]isoquinoline (N at 2; positions 1, 3-10 are CH)
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
    # benzo[h]isoquinoline (N at 6; positions 1-5, 7-10 are CH)
    ("c1ccc2c(c1)ccc1ccncc12",         "benzo[h]isoquinoline"),
    ("Cc1cc2ccncc2c2ccccc12",          "1-methylbenzo[h]isoquinoline"),
    ("Cc1cccc2c1ccc1ccncc12",          "2-methylbenzo[h]isoquinoline"),
    ("Cc1ccc2c(ccc3ccncc32)c1",        "3-methylbenzo[h]isoquinoline"),
    ("Cc1ccc2ccc3ccncc3c2c1",          "4-methylbenzo[h]isoquinoline"),
    ("Cc1cccc2ccc3ccncc3c12",          "5-methylbenzo[h]isoquinoline"),
    ("Cc1nccc2ccc3ccccc3c12",          "6-methylbenzo[h]isoquinoline"),
    ("Cc1cc2ccc3ccccc3c2cn1",          "8-methylbenzo[h]isoquinoline"),
    ("Cc1cncc2c1ccc1ccccc12",          "9-methylbenzo[h]isoquinoline"),
    ("Cc1cc2ccccc2c2cnccc12",          "10-methylbenzo[h]isoquinoline"),
])
def test_phase691(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
