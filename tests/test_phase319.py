"""Phase 319: cycloalkane dicarbothioamide (IUPAC 2013 P-65.1).

Extends Phase 306 dithioamide handler to support non-aromatic (cycloalkane) rings.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # cyclohexane dicarbothioamide
    ("NC(=S)C1CCC(C(=S)N)CC1",   "cyclohexane-1,4-dicarbothioamide"),
    ("NC(=S)C1CCCCC1C(=S)N",     "cyclohexane-1,2-dicarbothioamide"),
    # cyclopentane dicarbothioamide
    ("NC(=S)C1CCCC1C(=S)N",      "cyclopentane-1,2-dicarbothioamide"),
    # regressions: chain dithioamide unchanged
    ("NC(=S)C(=S)N",             "ethanedithioamide"),
    ("NC(=S)CC(=S)N",            "propanedithioamide"),
    # regression: benzene dicarbothioamide unchanged
    ("NC(=S)c1ccc(C(=S)N)cc1",  "benzene-1,4-dicarbothioamide"),
])
def test_phase319_cycloalkane_dithioamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
