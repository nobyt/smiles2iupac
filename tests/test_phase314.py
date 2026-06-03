"""Phase 314: cycloalkane carbothioamide/carboselenoamide (IUPAC 2013 P-66.3.1).

Exocyclic thioamide/selenoamide on a cycloalkane ring → "cyclo{stem}ane-carbothioamide".
Handlers now return None when adjacent to a ring, delegating to the cyclic path.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("NC(=S)C1CCCC1",   "cyclopentanecarbothioamide"),
    ("NC(=S)C1CCCCC1",  "cyclohexanecarbothioamide"),
    ("NC(=S)C1CCC1",    "cyclobutanecarbothioamide"),
    # regressions: open-chain thioamides unchanged
    ("CC(=S)N",         "ethanethioamide"),
    ("CCC(=S)N",        "propanethioamide"),
    ("NC(=S)N",         "urea-like or thiourea"),  # skip - not in scope
])
def test_phase314_cycloalkane_carbothioamide(smiles, expected):
    if expected == "urea-like or thiourea":
        pytest.skip("thiourea outside scope")
    assert smiles_to_iupac(smiles) == expected
