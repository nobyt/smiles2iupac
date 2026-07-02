"""Phase 771: benzo-naphthyridine α-ol/thiol → tautomers (IUPAC 2013).

- benzo[b][1,5]naphthyridine: C2 (alpha to N1) → 2(1H)-one
- benzo[b][1,6]naphthyridine: C1 and C3 (alpha to N2) → 1(2H)-one and 3(2H)-one
- benzo[b][1,8]naphthyridine: C2 (alpha to N1) → 2(1H)-one
- benzo[c][1,6]naphthyridine: C1,C3 (alpha to N2) and C6 (alpha to N5) → 1(2H)-one, 3(2H)-one, 6(5H)-one
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[b][1,5]naphthyridine C2
    ("Oc1ccc2nc3ccccc3cc2n1",       "benzo[b][1,5]naphthyridin-2(1H)-one"),
    ("Sc1ccc2nc3ccccc3cc2n1",       "benzo[b][1,5]naphthyridin-2(1H)-thione"),
    # benzo[b][1,6]naphthyridine C1 and C3
    ("Oc1nccc2nc3ccccc3cc12",       "benzo[b][1,6]naphthyridin-1(2H)-one"),
    ("Sc1nccc2nc3ccccc3cc12",       "benzo[b][1,6]naphthyridin-1(2H)-thione"),
    ("Oc1cc2nc3ccccc3cc2cn1",       "benzo[b][1,6]naphthyridin-3(2H)-one"),
    ("Sc1cc2nc3ccccc3cc2cn1",       "benzo[b][1,6]naphthyridin-3(2H)-thione"),
    # benzo[b][1,8]naphthyridine C2
    ("Oc1ccc2cc3ccccc3nc2n1",       "benzo[b][1,8]naphthyridin-2(1H)-one"),
    ("Sc1ccc2cc3ccccc3nc2n1",       "benzo[b][1,8]naphthyridin-2(1H)-thione"),
    # benzo[c][1,6]naphthyridine C6, C3, C1
    ("Oc1nc2ccncc2c2ccccc12",       "benzo[c][1,6]naphthyridin-6(5H)-one"),
    ("Sc1nc2ccncc2c2ccccc12",       "benzo[c][1,6]naphthyridin-6(5H)-thione"),
    ("Oc1cc2ncc3ccccc3c2cn1",       "benzo[c][1,6]naphthyridin-3(2H)-one"),
    ("Sc1cc2ncc3ccccc3c2cn1",       "benzo[c][1,6]naphthyridin-3(2H)-thione"),
    ("Oc1nccc2ncc3ccccc3c12",       "benzo[c][1,6]naphthyridin-1(2H)-one"),
    ("Sc1nccc2ncc3ccccc3c12",       "benzo[c][1,6]naphthyridin-1(2H)-thione"),
    # Regression: parent rings unaffected
    ("c1ccc2nc3ccccc3cc2n1",        "benzo[b][1,5]naphthyridine"),
    ("c1ccc2nc3ccncc3cc2c1",        "benzo[b][1,6]naphthyridine"),
    ("c1ccc2cc3ccccc3nc2n1",        "benzo[b][1,8]naphthyridine"),
    ("c1ccc2c(c1)cnc1ccncc12",      "benzo[c][1,6]naphthyridine"),
])
def test_phase771_benzo_naphthyridine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
