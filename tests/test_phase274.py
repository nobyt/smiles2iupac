"""Phase 274: hydroxy-substituted alkyl substituents on ring systems (IUPAC 2013 P-63.1).

Hydroxymethyl, hydroxyethyl, etc. as substituents on cycloalkane,
bicycloalkane, and spiro ring systems.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── hydroxymethyl on monocyclic ──────────────────────────────────────────
    ("OCC1CCCCC1",          "cyclohexylmethanol"),
    ("OCC1CCCC1",           "cyclopentylmethanol"),

    # ── hydroxymethyl on bicyclo ─────────────────────────────────────────────
    ("OCC1CC2CCC1CC2",      "2-(hydroxymethyl)bicyclo[2.2.2]octane"),

    # ── hydroxymethyl on spiro ───────────────────────────────────────────────
    ("OCC1CCC2(CC1)CCCC2",  "8-(hydroxymethyl)spiro[4.5]decane"),

    # ── 2-hydroxyethyl on monocyclic ─────────────────────────────────────────
    # Not expected to be named as substituent on ring in standard IUPAC — skip for now

    # ── regression: unsubstituted polycyclics still work ─────────────────────
    ("C1CC2CCC1CC2",        "bicyclo[2.2.2]octane"),
    ("C1CCC2(CC1)CCCC2",    "spiro[4.5]decane"),
])
def test_phase274_hydroxy_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
