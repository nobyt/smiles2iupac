"""Phase 757: 1,2,3-triazin-4-ol and 9H/7H-purin-2-ol → lactam tautomers (IUPAC 2013).

Extends triazine and purine tautomeric conversions:
- 1,2,3-triazin-4-ol → 1,2,3-triazin-4(3H)-one (and 4-thiol → 4(3H)-thione)
- 9H-purin-2-ol → 9H-purin-2(1H)-one
- 7H-purin-2-ol → 7H-purin-2(1H)-one
(and thiol → thione counterparts)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,2,3-triazin-4-ol → 1,2,3-triazin-4(3H)-one
    ("Oc1ccnnn1",                  "1,2,3-triazin-4(3H)-one"),
    ("Sc1ccnnn1",                  "1,2,3-triazin-4(3H)-thione"),
    # 9H-purin-2-ol → 9H-purin-2(1H)-one
    ("Oc1nc2[nH]cnc2cn1",         "9H-purin-2(1H)-one"),
    ("Oc1nc2nc[nH]c2cn1",         "7H-purin-2(1H)-one"),
    ("Sc1nc2[nH]cnc2cn1",         "9H-purin-2(1H)-thione"),
    ("Sc1nc2nc[nH]c2cn1",         "7H-purin-2(1H)-thione"),
    # Regression: parent rings unaffected
    ("c1ccnnn1",                  "1,2,3-triazine"),
    ("c1ncc2nc[nH]c2n1",         "9H-purine"),
    # Regression: Phase 755/756 purine-6 cases unchanged
    ("Oc1ncnc2[nH]cnc12",        "9H-purin-6(1H)-one"),
    ("Sc1ncnc2[nH]cnc12",        "9H-purin-6(1H)-thione"),
])
def test_phase757_triazine_and_purin2_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
