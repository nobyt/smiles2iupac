"""Phase 765: 2,7/2,6/1,6/3,6/4,5-phenanthroline α-ol/thiol → tautomers (IUPAC 2013).

Extends tautomeric conversions to remaining phenanthroline isomers:
- 2,7-phenanthroline: positions 1(N2), 6(N7), 8(N7)
- 2,6-phenanthroline: positions 1(N2), 3(N2), 5(N6)
- 1,6-phenanthroline: positions 2(N1), 5(N6)
- 3,6-phenanthroline: positions 2(N3), 4(N3), 5(N6)
- 4,5-phenanthroline: positions 3(N4), 6(N5)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2,7-phenanthroline
    ("Oc1ccc2c(ccc3ccncc32)n1",        "2,7-phenanthrolin-1(2H)-one"),
    ("Sc1ccc2c(ccc3ccncc32)n1",        "2,7-phenanthrolin-1(2H)-thione"),
    ("Oc1cc2ccc3ncccc3c2cn1",          "2,7-phenanthrolin-6(7H)-one"),
    ("Sc1cc2ccc3ncccc3c2cn1",          "2,7-phenanthrolin-6(7H)-thione"),
    ("Oc1nccc2ccc3ncccc3c12",          "2,7-phenanthrolin-8(7H)-one"),
    ("Sc1nccc2ccc3ncccc3c12",          "2,7-phenanthrolin-8(7H)-thione"),
    # 2,6-phenanthroline
    ("Oc1nccc2cnc3ccccc3c12",          "2,6-phenanthrolin-1(2H)-one"),
    ("Sc1nccc2cnc3ccccc3c12",          "2,6-phenanthrolin-1(2H)-thione"),
    ("Oc1cc2cnc3ccccc3c2cn1",          "2,6-phenanthrolin-3(2H)-one"),
    ("Sc1cc2cnc3ccccc3c2cn1",          "2,6-phenanthrolin-3(2H)-thione"),
    ("Oc1nc2ccccc2c2cnccc12",          "2,6-phenanthrolin-5(6H)-one"),
    ("Sc1nc2ccccc2c2cnccc12",          "2,6-phenanthrolin-5(6H)-thione"),
    # 1,6-phenanthroline
    ("Oc1ccc2cnc3ccccc3c2n1",          "1,6-phenanthrolin-2(1H)-one"),
    ("Sc1ccc2cnc3ccccc3c2n1",          "1,6-phenanthrolin-2(1H)-thione"),
    ("Oc1nc2ccccc2c2ncccc12",          "1,6-phenanthrolin-5(6H)-one"),
    ("Sc1nc2ccccc2c2ncccc12",          "1,6-phenanthrolin-5(6H)-thione"),
    # 3,6-phenanthroline
    ("Oc1cc2c(cn1)cnc1ccccc12",        "3,6-phenanthrolin-2(3H)-one"),
    ("Sc1cc2c(cn1)cnc1ccccc12",        "3,6-phenanthrolin-2(3H)-thione"),
    ("Oc1nccc2c1cnc1ccccc12",          "3,6-phenanthrolin-4(3H)-one"),
    ("Sc1nccc2c1cnc1ccccc12",          "3,6-phenanthrolin-4(3H)-thione"),
    ("Oc1nc2ccccc2c2ccncc12",          "3,6-phenanthrolin-5(6H)-one"),
    ("Sc1nc2ccccc2c2ccncc12",          "3,6-phenanthrolin-5(6H)-thione"),
    # 4,5-phenanthroline
    ("Oc1nc2ncccc2c2ccccc12",          "4,5-phenanthrolin-3(4H)-one"),
    ("Sc1nc2ncccc2c2ccccc12",          "4,5-phenanthrolin-3(4H)-thione"),
    ("Oc1ccc2c(ncc3ccccc32)n1",        "4,5-phenanthrolin-6(5H)-one"),
    ("Sc1ccc2c(ncc3ccccc32)n1",        "4,5-phenanthrolin-6(5H)-thione"),
    # Regression: parent rings unaffected
    ("c1cnc2ccc3ccncc3c2c1",           "2,7-phenanthroline"),
    ("c1ccc2c(c1)ncc1ccncc12",         "2,6-phenanthroline"),
    ("c1cnc2c(c1)cnc1ccccc12",         "1,6-phenanthroline"),
    ("c1ccc2c(c1)ncc1cnccc12",         "3,6-phenanthroline"),
    ("c1ccc2c(c1)cnc1ncccc12",         "4,5-phenanthroline"),
    # Regression: Phase 764 unchanged
    ("Oc1ccc2ccc3ncccc3c2n1",          "1,7-phenanthrolin-2(1H)-one"),
])
def test_phase765_phenanthroline_remaining_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
