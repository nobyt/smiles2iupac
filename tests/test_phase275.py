"""Phase 275: alkylidenecycloalkane — exo C=C double bond on ring (IUPAC 2013 P-31.1.6.1).

methylenecyclohexane, ethylidenecyclopentane, etc.
The exo =CH2 (and =CHR) group is named with the alkylidene prefix.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── methylenecycloalkane ───────────────────────────────────────────────
    ("C=C1CCCCC1",    "methylenecyclohexane"),
    ("C=C1CCCC1",     "methylenecyclopentane"),
    ("C=C1CCCCCC1",   "methylenecycloheptane"),

    # ── ethylidenecycloalkane ─────────────────────────────────────────────
    ("CC=C1CCCCC1",   "ethylidenecyclohexane"),
    ("CC=C1CCCC1",    "ethylidenecyclopentane"),

    # ── regressions: single bond methyl substituent still correct ─────────
    ("CC1CCCCC1",     "methylcyclohexane"),
    ("CC1CCCC1",      "methylcyclopentane"),
    # ketone exo =O unaffected
    ("O=C1CCCCC1",    "cyclohexanone"),
])
def test_phase275_alkylidenecycloalkane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
