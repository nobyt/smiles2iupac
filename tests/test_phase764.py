"""Phase 764: 1,7-phenanthroline and 1,8-phenanthroline α-ol/thiol → tautomers (IUPAC 2013).

Extends tautomeric conversions to asymmetric phenanthroline isomers:
- 1,7-phenanthrolin-2-ol → 1,7-phenanthrolin-2(1H)-one  (C2 alpha to N1)
- 1,7-phenanthrolin-8-ol → 1,7-phenanthrolin-8(7H)-one  (C8 alpha to N7)
- 1,8-phenanthrolin-2-ol → 1,8-phenanthrolin-2(1H)-one  (C2 alpha to N1)
- 1,8-phenanthrolin-9-ol → 1,8-phenanthrolin-9(8H)-one  (C9 alpha to N8)
- 1,8-phenanthrolin-7-ol → 1,8-phenanthrolin-7(8H)-one  (C7 alpha to N8)
(and thiol → thione counterparts)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,7-phenanthroline
    ("Oc1ccc2ccc3ncccc3c2n1",          "1,7-phenanthrolin-2(1H)-one"),
    ("Sc1ccc2ccc3ncccc3c2n1",          "1,7-phenanthrolin-2(1H)-thione"),
    ("Oc1ccc2c(ccc3cccnc32)n1",        "1,7-phenanthrolin-8(7H)-one"),
    ("Sc1ccc2c(ccc3cccnc32)n1",        "1,7-phenanthrolin-8(7H)-thione"),
    # 1,8-phenanthroline
    ("Oc1ccc2ccc3cnccc3c2n1",          "1,8-phenanthrolin-2(1H)-one"),
    ("Sc1ccc2ccc3cnccc3c2n1",          "1,8-phenanthrolin-2(1H)-thione"),
    ("Oc1cc2c(ccc3cccnc32)cn1",        "1,8-phenanthrolin-9(8H)-one"),
    ("Sc1cc2c(ccc3cccnc32)cn1",        "1,8-phenanthrolin-9(8H)-thione"),
    ("Oc1nccc2c1ccc1cccnc12",          "1,8-phenanthrolin-7(8H)-one"),
    ("Sc1nccc2c1ccc1cccnc12",          "1,8-phenanthrolin-7(8H)-thione"),
    # Regression: parent rings unaffected
    ("c1cnc2c(c1)ccc1ncccc12",         "1,7-phenanthroline"),
    ("c1cnc2c(c1)ccc1cnccc12",         "1,8-phenanthroline"),
    # Regression: Phase 763 unchanged
    ("Oc1ccc2ccc3cccnc3c2n1",          "1,10-phenanthrolin-2(1H)-one"),
])
def test_phase764_phenanthroline_asymmetric_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
