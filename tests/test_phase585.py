"""Phase 585: Substituted benzo[b][1,5]-, benzo[b][1,6]-, benzo[b][1,7]-,
benzo[b][1,8]-, and benzo[c][1,6]naphthyridine naming.
N positions vary by isomer; substitutable C positions: 8 per compound.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[b][1,7]naphthyridine (N at 2,10; sub C: 1,3-9)
    ("c1ccc2nc3cnccc3cc2c1",        "benzo[b][1,7]naphthyridine"),
    ("Cc1nccc2cc3ccccc3nc12",       "1-methylbenzo[b][1,7]naphthyridine"),
    ("Cc1cc2cc3ccccc3nc2cn1",       "3-methylbenzo[b][1,7]naphthyridine"),
    ("Cc1cncc2nc3ccccc3cc12",       "4-methylbenzo[b][1,7]naphthyridine"),
    ("Cc1c2ccccc2nc2cnccc12",       "5-methylbenzo[b][1,7]naphthyridine"),
    ("Cc1cccc2nc3cnccc3cc12",       "6-methylbenzo[b][1,7]naphthyridine"),
    ("Cc1ccc2nc3cnccc3cc2c1",       "7-methylbenzo[b][1,7]naphthyridine"),
    ("Cc1ccc2cc3ccncc3nc2c1",       "8-methylbenzo[b][1,7]naphthyridine"),
    ("Cc1cccc2cc3ccncc3nc12",       "9-methylbenzo[b][1,7]naphthyridine"),
    # benzo[b][1,8]naphthyridine (N at 1,10; sub C: 2-9)
    ("c1ccc2nc3ncccc3cc2c1",        "benzo[b][1,8]naphthyridine"),
    ("Cc1ccc2cc3ccccc3nc2n1",       "2-methylbenzo[b][1,8]naphthyridine"),
    ("Cc1cnc2nc3ccccc3cc2c1",       "3-methylbenzo[b][1,8]naphthyridine"),
    ("Cc1ccnc2nc3ccccc3cc12",       "4-methylbenzo[b][1,8]naphthyridine"),
    ("Cc1c2ccccc2nc2ncccc12",       "5-methylbenzo[b][1,8]naphthyridine"),
    ("Cc1cccc2nc3ncccc3cc12",       "6-methylbenzo[b][1,8]naphthyridine"),
    ("Cc1ccc2nc3ncccc3cc2c1",       "7-methylbenzo[b][1,8]naphthyridine"),
    ("Cc1ccc2cc3cccnc3nc2c1",       "8-methylbenzo[b][1,8]naphthyridine"),
    ("Cc1cccc2cc3cccnc3nc12",       "9-methylbenzo[b][1,8]naphthyridine"),
    # benzo[b][1,5]naphthyridine (N at 1,5; sub C: 2-4,6-10)
    ("c1ccc2nc3cccnc3cc2c1",        "benzo[b][1,5]naphthyridine"),
    ("Cc1ccc2nc3ccccc3cc2n1",       "2-methylbenzo[b][1,5]naphthyridine"),
    ("Cc1cnc2cc3ccccc3nc2c1",       "3-methylbenzo[b][1,5]naphthyridine"),
    ("Cc1ccnc2cc3ccccc3nc12",       "4-methylbenzo[b][1,5]naphthyridine"),
    ("Cc1cccc2cc3ncccc3nc12",       "6-methylbenzo[b][1,5]naphthyridine"),
    ("Cc1ccc2cc3ncccc3nc2c1",       "7-methylbenzo[b][1,5]naphthyridine"),
    ("Cc1ccc2nc3cccnc3cc2c1",       "8-methylbenzo[b][1,5]naphthyridine"),
    ("Cc1cccc2nc3cccnc3cc12",       "9-methylbenzo[b][1,5]naphthyridine"),
    ("Cc1c2ccccc2nc2cccnc12",       "10-methylbenzo[b][1,5]naphthyridine"),
    # benzo[b][1,6]naphthyridine (N at 2,5; sub C: 1,3,4,6-10)
    ("c1ccc2nc3ccncc3cc2c1",        "benzo[b][1,6]naphthyridine"),
    ("Cc1nccc2nc3ccccc3cc12",       "1-methylbenzo[b][1,6]naphthyridine"),
    ("Cc1cc2nc3ccccc3cc2cn1",       "3-methylbenzo[b][1,6]naphthyridine"),
    ("Cc1cncc2cc3ccccc3nc12",       "4-methylbenzo[b][1,6]naphthyridine"),
    ("Cc1cccc2cc3cnccc3nc12",       "6-methylbenzo[b][1,6]naphthyridine"),
    ("Cc1ccc2cc3cnccc3nc2c1",       "7-methylbenzo[b][1,6]naphthyridine"),
    ("Cc1ccc2nc3ccncc3cc2c1",       "8-methylbenzo[b][1,6]naphthyridine"),
    ("Cc1cccc2nc3ccncc3cc12",       "9-methylbenzo[b][1,6]naphthyridine"),
    ("Cc1c2ccccc2nc2ccncc12",       "10-methylbenzo[b][1,6]naphthyridine"),
    # benzo[c][1,6]naphthyridine (N at 2,5; sub C: 1,3,4,6-10)
    ("c1ccc2c(c1)cnc1ccncc12",      "benzo[c][1,6]naphthyridine"),
    ("Cc1nccc2ncc3ccccc3c12",       "1-methylbenzo[c][1,6]naphthyridine"),
    ("Cc1cc2ncc3ccccc3c2cn1",       "3-methylbenzo[c][1,6]naphthyridine"),
    ("Cc1cncc2c1ncc1ccccc12",       "4-methylbenzo[c][1,6]naphthyridine"),
    ("Cc1nc2ccncc2c2ccccc12",       "6-methylbenzo[c][1,6]naphthyridine"),
    ("Cc1cccc2c1cnc1ccncc12",       "7-methylbenzo[c][1,6]naphthyridine"),
    ("Cc1ccc2c(cnc3ccncc32)c1",     "8-methylbenzo[c][1,6]naphthyridine"),
    ("Cc1ccc2cnc3ccncc3c2c1",       "9-methylbenzo[c][1,6]naphthyridine"),
    ("Cc1cccc2cnc3ccncc3c12",       "10-methylbenzo[c][1,6]naphthyridine"),
])
def test_phase585_benzo_naphthyridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
