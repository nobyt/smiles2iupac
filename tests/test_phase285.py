"""Phase 285: thioketone naming — "propane-2-thione" (IUPAC 2013 P-31.1.2.2).

The suffix '-thione' begins with a consonant, so the terminal 'e' of the parent
hydride name is NOT elided: 'propane-2-thione', not 'propan-2-thione'.
(Contrast with '-one' which begins with 'o' and elides: 'propan-2-one'.)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CC(=S)C",    "propane-2-thione"),
    ("CCC(=S)C",   "butane-2-thione"),
    ("CCC(=S)CC",  "pentane-3-thione"),
    # ── regressions ───────────────────────────────────────────────────────
    ("CC=S",       "ethanethial"),
    ("CCC=S",      "propanethial"),
    ("S=C1CCCCC1", "cyclohexanethione"),
])
def test_phase285_thioketone_naming(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
