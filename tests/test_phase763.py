"""Phase 763: 1,10-phenanthroline and 4,7-phenanthroline α-ol/thiol → tautomers (IUPAC 2013).

Extends tautomeric conversions to symmetric phenanthroline isomers:
- 1,10-phenanthrolin-2-ol → 1,10-phenanthrolin-2(1H)-one
- 4,7-phenanthrolin-3-ol → 4,7-phenanthrolin-3(4H)-one
(and thiol → thione counterparts)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,10-phenanthroline: C2 alpha to N1
    ("Oc1ccc2ccc3cccnc3c2n1",          "1,10-phenanthrolin-2(1H)-one"),
    ("Sc1ccc2ccc3cccnc3c2n1",          "1,10-phenanthrolin-2(1H)-thione"),
    # 4,7-phenanthroline: C3 alpha to N4
    ("Oc1ccc2c(ccc3ncccc32)n1",        "4,7-phenanthrolin-3(4H)-one"),
    ("Sc1ccc2c(ccc3ncccc32)n1",        "4,7-phenanthrolin-3(4H)-thione"),
    # Regression: parent rings unaffected
    ("c1ccnc2c1ccc1cccnc12",           "1,10-phenanthroline"),
    ("c1cnc2ccc3ncccc3c2c1",           "4,7-phenanthroline"),
    # Regression: Phase 762 unchanged
    ("Oc1cccc2ccc3cccnc3c21",          "benzo[h]quinolin-1(10H)-one"),
])
def test_phase763_phenanthroline_symmetric_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
