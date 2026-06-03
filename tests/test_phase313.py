"""Phase 313: thiolactone naming (IUPAC 2013 P-65.6.3.4).

Cyclic thioesters (thiolactones) → "{ring}an-{loc}-one".
The thioester detector now skips when both S and C=O are in the same ring;
the heterocycle handler then applies the -one suffix.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("O=C1CCCS1",    "thiolan-2-one"),   # thiolane-2-one / thio-gamma-butyrolactone
    ("O=C1CCS1",     "thietan-2-one"),   # thiolanone 4-membered
    ("O=C1CCCCS1",   "thian-2-one"),     # thiane-2-one
    # regressions: open-chain thioesters unchanged
    ("CC(=O)SC",     "S-methyl ethanethioate"),
    ("CC(=O)SCC",    "S-ethyl ethanethioate"),
    # regressions: oxygen lactones unchanged
    ("O=C1CCCO1",    "oxolan-2-one"),
    ("O=C1CCCCO1",   "oxan-2-one"),
    # regressions: lactam unchanged
    ("O=C1CCCN1",    "pyrrolidin-2-one"),
])
def test_phase313_thiolactone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
