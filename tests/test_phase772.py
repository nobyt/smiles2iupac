"""Phase 772: 3,5-phenanthroline α-ol/thiol → tautomers (IUPAC 2013).

3,5-phenanthroline has N at positions 3 and 5.
- C2 (alpha to N3) → 2(3H)-one
- C4 (alpha to N5) → 4(5H)-one
- C6 (alpha to N5) → 6(5H)-one
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 3,5-phenanthroline hydroxy
    ("Oc1nc2cnccc2c2ccccc12",       "3,5-phenanthrolin-2(3H)-one"),
    ("Oc1nccc2c1ncc1ccccc12",       "3,5-phenanthrolin-4(5H)-one"),
    ("Oc1cc2c(cn1)ncc1ccccc12",     "3,5-phenanthrolin-6(5H)-one"),
    # 3,5-phenanthroline thiol
    ("Sc1nc2cnccc2c2ccccc12",       "3,5-phenanthrolin-2(3H)-thione"),
    ("Sc1nccc2c1ncc1ccccc12",       "3,5-phenanthrolin-4(5H)-thione"),
    ("Sc1cc2c(cn1)ncc1ccccc12",     "3,5-phenanthrolin-6(5H)-thione"),
    # Regression: parent ring unaffected
    ("c1ccc2c(c1)cnc1cnccc12",      "3,5-phenanthroline"),
    # Regression: Phase 765 unchanged
    ("Oc1nc2ncccc2c2ccccc12",       "4,5-phenanthrolin-3(4H)-one"),
])
def test_phase772_phenanthroline_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
