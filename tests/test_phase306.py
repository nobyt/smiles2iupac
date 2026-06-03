"""Phase 306: dithioamide naming (IUPAC 2013 P-65.1.2.4).

Two thioamide groups on a chain → "ethanedithioamide" etc.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("NC(=S)C(=S)N",       "ethanedithioamide"),
    ("NC(=S)CC(=S)N",      "propanedithioamide"),
    ("NC(=S)CCC(=S)N",     "butanedithioamide"),
    # regression: mono thioamide unchanged
    ("CC(=S)N",            "ethanethioamide"),
    ("C(=S)N",             "methanethioamide"),
])
def test_phase306_dithioamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
