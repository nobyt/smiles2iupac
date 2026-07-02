"""Phase 756: 7H-purin-6-ol, purine-6-thiol, 1,2,4-triazine-5/6-ol → tautomers (IUPAC 2013).

Extends purine and adds 1,2,4-triazine tautomeric conversions:
- 7H-purin-6-ol → 7H-purin-6(1H)-one
- 9H-purine-6-thiol → 9H-purin-6(1H)-thione  (6-mercaptopurine)
- 7H-purine-6-thiol → 7H-purin-6(1H)-thione
- 1,2,4-triazin-6-ol → 1,2,4-triazin-6(1H)-one
- 1,2,4-triazin-5-ol → 1,2,4-triazin-5(4H)-one
(and thione counterparts)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 7H-purin-6-ol → 7H-purin-6(1H)-one
    ("Oc1ncnc2nc[nH]c12",          "7H-purin-6(1H)-one"),
    # purine-6-thiol tautomers (6-mercaptopurine)
    ("Sc1ncnc2[nH]cnc12",          "9H-purin-6(1H)-thione"),
    ("Sc1ncnc2nc[nH]c12",          "7H-purin-6(1H)-thione"),
    # 1,2,4-triazine tautomers
    ("Oc1cncnn1",                  "1,2,4-triazin-6(1H)-one"),
    ("Oc1cnncn1",                  "1,2,4-triazin-5(4H)-one"),
    ("Sc1cncnn1",                  "1,2,4-triazin-6(1H)-thione"),
    ("Sc1cnncn1",                  "1,2,4-triazin-5(4H)-thione"),
    # Regression: Phase 755 purine-6-ol unchanged
    ("Oc1ncnc2[nH]cnc12",          "9H-purin-6(1H)-one"),
    # Regression: parent rings unaffected
    ("c1ncc2nc[nH]c2n1",           "9H-purine"),
    ("c1ncc2[nH]cnc2n1",           "7H-purine"),
    ("c1cncnn1",                   "1,2,4-triazine"),
])
def test_phase756_purine_thiol_and_triazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
