"""Phase 284: imine locant always cited for chains ≥ 2 carbons (IUPAC 2013 P-62.3.1.1).

IUPAC 2013 requires the locant before "-imine" for all chains with ≥ 2 carbons,
matching the convention used for "-amine" and "-ol".  Only methanimine (1C) is
cited without a locant.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── 1C: no locant (unambiguous) ───────────────────────────────────────
    ("C=N",            "methanimine"),
    ("C=NC",           "N-methylmethanimine"),

    # ── 2C: locant required ───────────────────────────────────────────────
    ("CC=N",           "ethanimine"),
    ("CC=NC",          "N-methylethanimine"),
    ("CC=NCC",         "N-ethylethanimine"),

    # ── 3C+: locant required ─────────────────────────────────────────────
    ("CCC=N",          "propan-1-imine"),
    ("CC(=N)C",        "propan-2-imine"),
    ("CCC=NC",         "N-methylpropan-1-imine"),

    # ── regressions: other suffixes unaffected ────────────────────────────
    ("CC=NN",          "ethanal hydrazone"),
    ("CC=NO",          "ethanal oxime"),
    ("C=CC=N",         "prop-2-en-1-imine"),
])
def test_phase284_imine_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
