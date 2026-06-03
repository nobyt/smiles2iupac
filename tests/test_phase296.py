"""Phase 296: selenoamide naming (IUPAC 2013 P-66.8.3 / P-65.1.1.4).

C(=[Se])-NH2 → {alkane}selenoamide  (analog of thioamide with Se)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # basic selenoamides
    ("[Se]=CN",            "methaneselenoamide"),
    ("CC(=[Se])N",         "ethaneselenoamide"),
    ("CCC(=[Se])N",        "propaneselenoamide"),
    # N-substituted
    ("CC(=[Se])NC",        "N-methylethaneselenoamide"),
    # regressions: thioamide unchanged
    ("CC(=S)N",            "ethanethioamide"),
    ("C(=S)N",             "methanethioamide"),
])
def test_phase296_selenoamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
