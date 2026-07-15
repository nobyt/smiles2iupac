"""Phase 298: telluramide naming (IUPAC 2013 P-66.8.3 / P-65.1.1.4).

C(=[Te])-NH2 → {alkane}telluramide  (analog of thioamide / selenoamide with Te)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("[Te]=CN",            "methanetelluramide"),
    ("CC(=[Te])N",         "ethanetelluramide"),
    ("CCC(=[Te])N",        "propanetelluramide"),
    ("CC(=[Te])NC",        "N-methylethanetelluramide"),
    # regressions
    ("CC(=[Se])N",         "ethaneselenoamide"),
    ("CC(=S)N",            "ethanethioamide"),
])
def test_phase298_telluramide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
