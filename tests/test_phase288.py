"""Phase 288: polycyclic carboxylic acid suffix — "-carboxylic acid" not "carboxy-" prefix
(IUPAC 2013 P-65.1.1.1).

A carboxyl group on a polycyclic ring is expressed as "-carboxylic acid" suffix on
the ring name, with the terminal "e" of "ane"/"ene" retained (e.g., "octane-2-").
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── saturated bicyclic ────────────────────────────────────────────────
    ("OC(=O)C1CC2CCC1CC2",    "bicyclo[2.2.2]octane-2-carboxylic acid"),
    ("OC(=O)C1CC2CC1CC2",     "bicyclo[2.2.1]heptane-2-carboxylic acid"),

    # ── unsaturated bicyclic ──────────────────────────────────────────────
    ("OC(=O)C1CC2CC1C=C2",    "bicyclo[2.2.1]hept-5-ene-2-carboxylic acid"),

    # ── regressions ───────────────────────────────────────────────────────
    ("OC1CC2CCC1CC2",          "bicyclo[2.2.2]octan-2-ol"),
    ("C1CC2CCC1CC2",           "bicyclo[2.2.2]octane"),
    ("OC(=O)C1CCCCC1",         "cyclohexanecarboxylic acid"),
])
def test_phase288_polycyclic_carboxylic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
