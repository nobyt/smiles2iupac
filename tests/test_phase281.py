"""Phase 281: isocyanide substitutive PIN — "isocyano" prefix (IUPAC 2013 P-62.5.3.2).

The preferred IUPAC 2013 name for isocyanides uses the prefix "isocyano" on the
parent hydrocarbon.  The functional-class retained name "alkyl isocyanide" is
acceptable but is not the PIN.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── acyclic isocyanides ───────────────────────────────────────────────
    ("C[N+]#[C-]",          "isocyanomethane"),
    ("CC[N+]#[C-]",         "isocyanoethane"),
    ("[C-]#[N+]CCC",        "isocyanopropane"),
    ("[C-]#[N+]CCCC",       "isocyanobutane"),
    # ── aromatic isocyanide ───────────────────────────────────────────────
    ("c1ccccc1[N+]#[C-]",   "isocyanobenzene"),

    # ── regressions: nitrile and isocyanate unchanged ─────────────────────
    ("CC#N",                "acetonitrile"),
    ("CCC#N",               "propanenitrile"),
    ("CN=C=O",              "isocyanatomethane"),
    ("CCN=C=O",             "isocyanatoethane"),
])
def test_phase281_isocyanide_pin(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
