"""Phase 775: missing naphthyridine/cinnoline α positions (IUPAC 2013).

1,6-naphthyridine (N1, N6): C5 and C7 are alpha to N6.
1,7-naphthyridine (N1, N7): C6 and C8 are alpha to N7.
cinnoline (N1, N2): C3 is alpha to N2 (C4 was already handled by Phase 748).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,6-naphthyridine: C5 and C7 alpha to N6
    ("Oc1nccc2ncccc12",     "1,6-naphthyridin-5(6H)-one"),
    ("Sc1nccc2ncccc12",     "1,6-naphthyridin-5(6H)-thione"),
    ("Oc1cc2ncccc2cn1",     "1,6-naphthyridin-7(6H)-one"),
    ("Sc1cc2ncccc2cn1",     "1,6-naphthyridin-7(6H)-thione"),
    # 1,7-naphthyridine: C6 and C8 alpha to N7
    ("Oc1cc2cccnc2cn1",     "1,7-naphthyridin-6(7H)-one"),
    ("Sc1cc2cccnc2cn1",     "1,7-naphthyridin-6(7H)-thione"),
    ("Oc1nccc2cccnc12",     "1,7-naphthyridin-8(7H)-one"),
    ("Sc1nccc2cccnc12",     "1,7-naphthyridin-8(7H)-thione"),
    # cinnoline: C3 alpha to N2
    ("Oc1cc2ccccc2nn1",     "cinnolin-3(2H)-one"),
    ("Sc1cc2ccccc2nn1",     "cinnolin-3(2H)-thione"),
    # Regression: previously handled positions unchanged
    ("Oc1ccc2cnccc2n1",     "1,6-naphthyridin-2(1H)-one"),
    ("Oc1ccc2ccncc2n1",     "1,7-naphthyridin-2(1H)-one"),
    ("Oc1cnnc2ccccc12",     "cinnolin-4(1H)-one"),
    ("c1cc2ccccc2nn1",      "cinnoline"),
])
def test_phase775_missing_alpha_positions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
