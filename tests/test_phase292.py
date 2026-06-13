"""Phase 292: isocyanate/isothiocyanate functional-class PIN (IUPAC 2013 P-65.3.1).

Alkyl isocyanate/isothiocyanate functional class names are the PINs:
  CN=C=O  → methyl isocyanate
  CN=C=S  → methyl isothiocyanate
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── isocyanate PIN ────────────────────────────────────────────────────
    ("CN=C=O",          "methyl isocyanate"),
    ("CCN=C=O",         "ethyl isocyanate"),
    ("CCCN=C=O",        "propyl isocyanate"),
    ("CC(C)N=C=O",      "propan-2-yl isocyanate"),

    # ── isothiocyanate PIN ────────────────────────────────────────────────
    ("CN=C=S",          "methyl isothiocyanate"),
    ("CCN=C=S",         "ethyl isothiocyanate"),
    ("CCCN=C=S",        "propyl isothiocyanate"),

    # ── aromatic ──────────────────────────────────────────────────────────
    ("O=C=Nc1ccccc1",   "phenyl isocyanate"),
    ("S=C=Nc1ccccc1",   "phenyl isothiocyanate"),

    # ── regressions ──────────────────────────────────────────────────────
    ("CC#N",            "acetonitrile"),
])
def test_phase292_isocyanate_pin(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
