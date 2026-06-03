"""Phase 283: sulfoxide, sulfone, and disulfide substitutive PINs (IUPAC 2013 P-63.6–P-63.7).

Sulfoxides (R-SO-R'), sulfones (R-SO2-R'), and disulfides (R-SS-R') are named by
substitutive nomenclature: the shorter chain forms the sulfinyl/sulfonyl/disulfanyl
prefix while the longer chain is the parent.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── sulfoxides ────────────────────────────────────────────────────────
    ("CS(=O)C",              "dimethyl sulfoxide"),
    ("CCS(=O)CC",            "diethyl sulfoxide"),
    ("CCS(=O)C",             "ethyl methyl sulfoxide"),

    # ── sulfones ──────────────────────────────────────────────────────────
    ("CS(=O)(=O)C",          "dimethyl sulfone"),
    ("CCS(=O)(=O)CC",        "diethyl sulfone"),
    ("CCS(=O)(=O)C",         "ethyl methyl sulfone"),

    # ── disulfides ────────────────────────────────────────────────────────
    ("CSSC",                 "dimethyl disulfide"),
    ("CCSSCC",               "diethyl disulfide"),

    # ── regressions: sulfonamide, sulfonic acid, thiol unchanged ──────────
    ("CS(=O)(=O)N",          "methanesulfonamide"),
    ("CS(=O)(=O)O",          "methanesulfonic acid"),
    ("CCS",                  "ethanethiol"),
])
def test_phase283_sulfoxide_sulfone_disulfide_pin(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
