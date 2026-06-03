"""Phase 285: thioketone locant elision — "propan-2-thione" not "propane-2-thione"
(IUPAC 2013 P-31.1.2.2).

The terminal 'e' of a parent hydride name is elided before a locant-bearing suffix
(propane → propan- before -2-thione), matching the convention used for "-one".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CC(=S)C",    "propan-2-thione"),
    ("CCC(=S)C",   "butan-2-thione"),
    ("CCC(=S)CC",  "pentan-3-thione"),
    # ── regressions ───────────────────────────────────────────────────────
    ("CC=S",       "ethanethial"),
    ("CCC=S",      "propanethial"),
    ("S=C1CCCCC1", "cyclohexanethione"),
])
def test_phase285_thioketone_locant_elision(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
