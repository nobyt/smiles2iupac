"""Phase 286: selenide, diselenide, telluride substitutive PINs (IUPAC 2013 P-63.6).

Same substitutive pattern as sulfide (Phase 282): the shorter chain forms the
(R)selanyl / (R)diselanyl / (R)tellanyl prefix; longer chain is the parent.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── selenides ──────────────────────────────────────────────────────────
    ("C[Se]C",           "(methylselanyl)methane"),
    ("CC[Se]C",          "(methylselanyl)ethane"),
    ("CC[Se]CC",         "(ethylselanyl)ethane"),

    # ── diselenides ────────────────────────────────────────────────────────
    ("C[Se][Se]C",       "(methyldiselanyl)methane"),
    ("CC[Se][Se]CC",     "(ethyldiselanyl)ethane"),

    # ── tellurides ────────────────────────────────────────────────────────
    ("C[Te]C",           "(methyltellanyl)methane"),

    # ── regressions ───────────────────────────────────────────────────────
    ("C[SeH]",           "methaneselenol"),
    ("CSC",              "(methylsulfanyl)methane"),
])
def test_phase286_selenide_telluride_pin(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
