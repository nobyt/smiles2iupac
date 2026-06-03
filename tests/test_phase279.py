"""Phase 279: adamantane and cubane retained cage compound names (IUPAC 2013 P-31.1.6.4).

Adamantane (tricyclo[3.3.1.1^{3,7}]decane) and cubane are retained names in IUPAC 2013.
Substituted adamantane: position 1 = bridgehead (CH), position 2 = bridge (CH2).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── unsubstituted cage compounds ────────────────────────────────────
    ("C1C2CC3CC1CC(C2)C3",        "adamantane"),
    ("C12C3C4C1C5C4C3C25",        "cubane"),

    # ── substituted adamantane ─────────────────────────────────────────
    # methyl at bridgehead (position 1)
    ("C1C2(C)CC3CC1CC(C2)C3",     "1-methyladamantane"),
    # methyl at bridge (position 2)
    ("CC1C2CC3CC1CC(C2)C3",       "2-methyladamantane"),

    # ── regressions: bicyclic and spiro still work ────────────────────
    ("C1CC2CCC1CC2",              "bicyclo[2.2.2]octane"),
    ("C1CC2CCC1C2",               "bicyclo[2.2.1]heptane"),
    ("C1CCC2(CC1)CCCC2",          "spiro[4.5]decane"),
])
def test_phase279_cage_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
