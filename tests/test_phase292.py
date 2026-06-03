"""Phase 292: isocyanate/isothiocyanate substitutive PIN (IUPAC 2013 P-65.3.1).

Alkyl isocyanate/isothiocyanate functional class names are retained acceptable
but not PINs.  The PIN uses the substitutive prefix form:
  CN=C=O  → isocyanatomethane       (not methyl isocyanate)
  CN=C=S  → isothiocyanatomethane   (not methyl isothiocyanate)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── isocyanate PIN ────────────────────────────────────────────────────
    ("CN=C=O",          "isocyanatomethane"),
    ("CCN=C=O",         "isocyanatoethane"),
    ("CCCN=C=O",        "isocyanatopropane"),
    ("CC(C)N=C=O",      "2-isocyanatopropane"),

    # ── isothiocyanate PIN ────────────────────────────────────────────────
    ("CN=C=S",          "isothiocyanatomethane"),
    ("CCN=C=S",         "isothiocyanatoethane"),
    ("CCCN=C=S",        "isothiocyanatopropane"),

    # ── aromatic: already substitutive, regression ─────────────────────────
    ("O=C=Nc1ccccc1",   "isocyanatobenzene"),
    ("S=C=Nc1ccccc1",   "isothiocyanatobenzene"),

    # ── regressions ──────────────────────────────────────────────────────
    ("CC#N",            "acetonitrile"),
])
def test_phase292_isocyanate_pin(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
