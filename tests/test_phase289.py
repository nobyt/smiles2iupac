"""Phase 289: unsaturated polycyclic + PCG — PCG gets lower locant than double bond
(IUPAC 2013 P-44.1.1).

When a bicyclic ring has both a double bond and a principal characteristic group
(PCG), the PCG locant is minimized first, then double bond locants.  The PCG
is expressed as a suffix on the ring name.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── norbornene-ol ─────────────────────────────────────────────────────
    ("OC1CC2CC1C=C2",   "bicyclo[2.2.1]hept-5-en-2-ol"),

    # ── norbornenone ──────────────────────────────────────────────────────
    # C1C2CC(=O)C1C=C2 — ring ketone + double bond
    ("O=C1CC2CC1C=C2",  "bicyclo[2.2.1]hept-5-en-2-one"),

    # ── regressions: no PCG → double bond gets lowest locant ──────────────
    ("C1CC2CC1C=C2",    "bicyclo[2.2.1]hept-2-ene"),
    ("C1C2C=CC1C=C2",   "bicyclo[2.2.1]hepta-2,5-diene"),
])
def test_phase289_unsaturated_bicyclic_pcg(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
