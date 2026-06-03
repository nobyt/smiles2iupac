"""Phase 287: polycyclic compounds with principal characteristic groups use suffix
nomenclature (IUPAC 2013 P-31.1.3, P-44.1).

A single PCG (hydroxy → -ol, oxo → -one, amino → -amine) on a saturated
bicyclic or spiro ring is expressed as a suffix on the ring name rather than
as a hydroxy/oxo/amino prefix.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── bicyclo -ol ────────────────────────────────────────────────────────
    ("OC1CC2CCC1CC2",   "bicyclo[2.2.2]octan-2-ol"),
    ("OC1CC2CC1CC2",    "bicyclo[2.2.1]heptan-2-ol"),

    # ── bicyclo -one ───────────────────────────────────────────────────────
    ("O=C1CC2CCC1CC2",  "bicyclo[2.2.2]octan-2-one"),
    ("O=C1CC2CC1CC2",   "bicyclo[2.2.1]heptan-2-one"),

    # ── bicyclo -amine ─────────────────────────────────────────────────────
    ("NC1CC2CCC1CC2",   "bicyclo[2.2.2]octan-2-amine"),

    # ── spiro -ol ──────────────────────────────────────────────────────────
    ("OC1CCC2(CC1)CCCC2",  "spiro[4.5]decan-8-ol"),

    # ── regressions: unsubstituted bicyclics unchanged ─────────────────────
    ("C1CC2CCC1CC2",    "bicyclo[2.2.2]octane"),
    ("C1CC2CC1CC2",     "bicyclo[2.2.1]heptane"),
    ("C1CCC2(CC1)CCCC2","spiro[4.5]decane"),
])
def test_phase287_polycyclic_pcg_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
