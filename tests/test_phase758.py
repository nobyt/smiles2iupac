"""Phase 758: quinazoline-2,4-diol and quinoxaline-2,3-diol → di-lactam tautomers (IUPAC 2013).

Extends the pyrimidine-2,4-diol → dione pattern (Phase 752) to benzo-diazines:
- quinazoline-2,4-diol → quinazoline-2,4(1H,3H)-dione
- quinoxaline-2,3-diol → quinoxaline-2,3(1H,4H)-dione
(and dithiol → dithione counterparts)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # quinazoline-2,4-diol → quinazoline-2,4(1H,3H)-dione
    ("Oc1nc(O)c2ccccc2n1",          "quinazoline-2,4(1H,3H)-dione"),
    ("Sc1nc(S)c2ccccc2n1",          "quinazoline-2,4(1H,3H)-dithione"),
    # quinoxaline-2,3-diol → quinoxaline-2,3(1H,4H)-dione
    ("Oc1nc2ccccc2nc1O",            "quinoxaline-2,3(1H,4H)-dione"),
    ("Sc1nc2ccccc2nc1S",            "quinoxaline-2,3(1H,4H)-dithione"),
    # Regression: mono-OH unchanged (Phase 743/748)
    ("Oc1nc2ccccc2nc1",             "quinoxalin-2(1H)-one"),
    ("Oc1ncnc2ccccc12",             "quinazolin-4(3H)-one"),
    # Regression: pyrimidine-2,4-diol unchanged (Phase 752)
    ("Oc1ccnc(O)n1",                "pyrimidine-2,4(1H,3H)-dione"),
    # Regression: parent rings unaffected
    ("c1ncnc2ccccc21",              "quinazoline"),
    ("c1cnc2ccccc2n1",              "quinoxaline"),
])
def test_phase758_quinazoline_quinoxaline_diol_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
