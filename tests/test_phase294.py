"""Phase 294: disulfonic acid / disulfonamide (IUPAC 2013 P-65.3.1).

Two sulfonic acid or sulfonamide groups on a chain: aggregate into "disulfonic acid"
/ "disulfonamide" with locants.  The suffix starts with 'd' (consonant) → terminal
'e' of parent is retained: "ethane-1,2-disulfonic acid" (not "ethan-1,2-...").
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── disulfonic acid ───────────────────────────────────────────────────
    ("OS(=O)(=O)CCS(=O)(=O)O",   "ethane-1,2-disulfonic acid"),
    ("OS(=O)(=O)CCCS(=O)(=O)O",  "propane-1,3-disulfonic acid"),

    # ── disulfonamide ─────────────────────────────────────────────────────
    ("NS(=O)(=O)CCS(=O)(=O)N",   "ethane-1,2-disulfonamide"),
    ("NS(=O)(=O)CCCS(=O)(=O)N",  "propane-1,3-disulfonamide"),

    # ── regressions: mono forms unchanged ────────────────────────────────
    ("NS(=O)(=O)C",              "methanesulfonamide"),
    ("OS(=O)(=O)C",              "methanesulfonic acid"),
    ("NS(=O)(=O)CC",             "ethanesulfonamide"),
    ("OS(=O)(=O)CC",             "ethanesulfonic acid"),
])
def test_phase294_disulfonic_disulfonamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
