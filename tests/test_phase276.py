"""Phase 276: unsaturated bicyclo[l.m.n]alkene/diene (IUPAC 2013 P-31.1.6.2).

Bicyclic ring systems with one or two double bonds within the ring.
Locant minimization: double bond locants take priority over substituent locants.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── bicyclo[2.2.1]hept-2-ene (norbornene) ────────────────────────────
    ("C1CC2CC1C=C2",           "bicyclo[2.2.1]hept-2-ene"),

    # ── bicyclo[2.2.2]oct-2-ene ──────────────────────────────────────────
    ("C1=CC2CCC1CC2",          "bicyclo[2.2.2]oct-2-ene"),

    # ── bicyclo[2.2.1]hepta-2,5-diene (norbornadiene) ────────────────────
    ("C1C2C=CC1C=C2",          "bicyclo[2.2.1]hepta-2,5-diene"),

    # ── substituted unsaturated bicyclic ─────────────────────────────────
    # 5-methyl on norbornene: methyl on the 1-carbon bridge
    ("CC1CC2CC1C=C2",          "5-methylbicyclo[2.2.1]hept-2-ene"),

    # ── regressions: saturated bicyclics unaffected ───────────────────────
    ("C1CC2CCC1CC2",           "bicyclo[2.2.2]octane"),
    ("C1CC2CCC1C2",            "bicyclo[2.2.1]heptane"),
    ("CC1CC2CCC1CC2",          "2-methylbicyclo[2.2.2]octane"),
])
def test_phase276_unsaturated_bicyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
