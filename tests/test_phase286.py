"""Phase 286: selenide, diselenide, telluride substitutive PINs (IUPAC 2013 P-63.6).

Same substitutive pattern as sulfide (Phase 282): the shorter chain forms the
(R)selanyl / (R)diselanyl / (R)tellanyl prefix; longer chain is the parent.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── selenides ──────────────────────────────────────────────────────────
    ("C[Se]C",           "dimethyl selenide"),
    ("CC[Se]C",          "ethyl methyl selenide"),
    ("CC[Se]CC",         "diethyl selenide"),

    # ── diselenides ────────────────────────────────────────────────────────
    ("C[Se][Se]C",       "dimethyl diselenide"),
    ("CC[Se][Se]CC",     "diethyl diselenide"),

    # ── tellurides ────────────────────────────────────────────────────────
    ("C[Te]C",           "dimethyl telluride"),

    # ── regressions ───────────────────────────────────────────────────────
    ("C[SeH]",           "methaneselenol"),
    ("CSC",              "dimethyl sulfide"),
])
def test_phase286_selenide_telluride_pin(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
