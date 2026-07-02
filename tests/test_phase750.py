"""Phase 750: hydroxy naphthyridines → preferred lactam tautomers (IUPAC 2013).

α-hydroxy positions (2 and 4 for 1,x-naphthyridines; 1 and 3 for 2,x-) adjacent
to ring N prefer the N-H keto (lactam) tautomer with indicated H at the proximal N.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,x-naphthyridine: 2-OH → 2(1H)-one
    ("Oc1nc2cccnc2cc1",          "1,5-naphthyridin-2(1H)-one"),
    ("Oc1nc2ccncc2cc1",          "1,6-naphthyridin-2(1H)-one"),
    ("Oc1nc2cnccc2cc1",          "1,7-naphthyridin-2(1H)-one"),
    ("Oc1nc2ncccc2cc1",          "1,8-naphthyridin-2(1H)-one"),
    # 1,x-naphthyridine: 4-OH → 4(1H)-one
    ("Oc1ccnc2cccnc12",          "1,5-naphthyridin-4(1H)-one"),
    ("Oc1ccnc2ccncc12",          "1,6-naphthyridin-4(1H)-one"),
    ("Oc1ccnc2cnccc12",          "1,7-naphthyridin-4(1H)-one"),
    ("Oc1ccnc2ncccc12",          "1,8-naphthyridin-4(1H)-one"),
    # Regression: parent naphthyridines unaffected
    ("c1cnc2cccnc2c1",           "1,5-naphthyridine"),
    ("c1cnc2ccncc2c1",           "1,6-naphthyridine"),
    ("c1cnc2cnccc2c1",           "1,7-naphthyridine"),
    ("c1cnc2ncccc2c1",           "1,8-naphthyridine"),
    ("c1cc2cnccc2cn1",           "2,6-naphthyridine"),
    ("c1cc2ccncc2cn1",           "2,7-naphthyridine"),
    # Regression: Phase 745 quinoline unchanged
    ("Oc1nc2ccccc2cc1",          "quinolin-2(1H)-one"),
    ("Oc1ccnc2ccccc12",          "quinolin-4(1H)-one"),
])
def test_phase750_naphthyridine_ol_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
