"""Phase 766: benzo[h/f/g]isoquinoline, benzo[g]quinoxaline, benzo[b][1,7]naphthyridine
α-ol/thiol → tautomers (IUPAC 2013).

- benzo[h]isoquinoline: positions 6(N7), 8(N7)
- benzo[f]isoquinoline: positions 5(N6), 7(N6)
- benzo[g]isoquinoline: positions 1(N2), 3(N2)
- benzo[g]quinoxaline: position 6(N5)
- benzo[b][1,7]naphthyridine: positions 1(N2), 3(N2)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[h]isoquinoline
    ("Oc1nccc2ccc3ccccc3c12",          "benzo[h]isoquinolin-6(7H)-one"),
    ("Sc1nccc2ccc3ccccc3c12",          "benzo[h]isoquinolin-6(7H)-thione"),
    ("Oc1cc2ccc3ccccc3c2cn1",          "benzo[h]isoquinolin-8(7H)-one"),
    ("Sc1cc2ccc3ccccc3c2cn1",          "benzo[h]isoquinolin-8(7H)-thione"),
    # benzo[f]isoquinoline
    ("Oc1cc2c(ccc3ccccc32)cn1",        "benzo[f]isoquinolin-5(6H)-one"),
    ("Sc1cc2c(ccc3ccccc32)cn1",        "benzo[f]isoquinolin-5(6H)-thione"),
    ("Oc1nccc2c1ccc1ccccc12",          "benzo[f]isoquinolin-7(6H)-one"),
    ("Sc1nccc2c1ccc1ccccc12",          "benzo[f]isoquinolin-7(6H)-thione"),
    # benzo[g]isoquinoline
    ("Oc1nccc2cc3ccccc3cc12",          "benzo[g]isoquinolin-1(2H)-one"),
    ("Sc1nccc2cc3ccccc3cc12",          "benzo[g]isoquinolin-1(2H)-thione"),
    ("Oc1cc2cc3ccccc3cc2cn1",          "benzo[g]isoquinolin-3(2H)-one"),
    ("Sc1cc2cc3ccccc3cc2cn1",          "benzo[g]isoquinolin-3(2H)-thione"),
    # benzo[g]quinoxaline
    ("Oc1cnc2cc3ccccc3cc2n1",          "benzo[g]quinoxalin-6(5H)-one"),
    ("Sc1cnc2cc3ccccc3cc2n1",          "benzo[g]quinoxalin-6(5H)-thione"),
    # benzo[b][1,7]naphthyridine
    ("Oc1nccc2cc3ccccc3nc12",          "benzo[b][1,7]naphthyridin-1(2H)-one"),
    ("Sc1nccc2cc3ccccc3nc12",          "benzo[b][1,7]naphthyridin-1(2H)-thione"),
    ("Oc1cc2cc3ccccc3nc2cn1",          "benzo[b][1,7]naphthyridin-3(2H)-one"),
    ("Sc1cc2cc3ccccc3nc2cn1",          "benzo[b][1,7]naphthyridin-3(2H)-thione"),
    # Regression: parent rings unaffected
    ("c1ccc2c(c1)ccc1ccncc12",         "benzo[h]isoquinoline"),
    ("c1ccc2c(c1)ccc1cnccc12",         "benzo[f]isoquinoline"),
    ("c1ccc2cc3cnccc3cc2c1",           "benzo[g]isoquinoline"),
    ("c1ccc2cc3nccnc3cc2c1",           "benzo[g]quinoxaline"),
    ("c1ccc2nc3cnccc3cc2c1",           "benzo[b][1,7]naphthyridine"),
    # Regression: Phase 765 unchanged
    ("Oc1ccc2c(ccc3ccncc32)n1",        "2,7-phenanthrolin-1(2H)-one"),
])
def test_phase766_benzo_isoquinoline_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
