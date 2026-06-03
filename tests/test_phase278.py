"""Phase 278: heteroatom bicyclo[l.m.n]alkane — "a" nomenclature (IUPAC 2013 P-31.1.2.3).

Bridged bicyclics with N, O, or S replacing ring carbons.
Heteroatom locants are minimized first, before double bond and substituent locants.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── DABCO: both bridgeheads are N ────────────────────────────────────
    ("C1CN2CCN1CC2",   "1,4-diazabicyclo[2.2.2]octane"),

    # ── ABCO: one bridgehead is N ─────────────────────────────────────────
    ("C1CN2CCC1CC2",   "1-azabicyclo[2.2.2]octane"),

    # ── regressions: pure carbocyclic still works ─────────────────────────
    ("C1CC2CCC1CC2",   "bicyclo[2.2.2]octane"),
    ("C1CC2CCC1C2",    "bicyclo[2.2.1]heptane"),
    ("C1CC2CC1C=C2",   "bicyclo[2.2.1]hept-2-ene"),
    ("C1C2C=CC1C=C2",  "bicyclo[2.2.1]hepta-2,5-diene"),
])
def test_phase278_hetero_bicyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
