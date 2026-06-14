"""Phase 293: cyanate/thiocyanate/azide substitutive PIN (IUPAC 2013 P-65.3.1).

Functional class names "methyl cyanate", "methyl thiocyanate", "methyl azide"
are retained acceptable but not PINs.  The PIN uses the substitutive prefix:
  COC#N        → cyanatomethane        (not methyl cyanate)
  CSC#N        → thiocyanatomethane    (not methyl thiocyanate)
  CN=[N+]=[N-] → azidomethane         (not methyl azide)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── cyanate PIN ───────────────────────────────────────────────────────
    ("COC#N",          "cyanatomethane"),
    ("CCOC#N",         "cyanatoethane"),
    ("CCCOC#N",        "cyanatopropane"),

    # ── thiocyanate PIN ──────────────────────────────────────────────────
    ("CSC#N",          "thiocyanatomethane"),
    ("CCSC#N",         "thiocyanatoethane"),
    ("CCCSC#N",        "thiocyanatopropane"),

    # ── azide PIN ────────────────────────────────────────────────────────
    ("CN=[N+]=[N-]",   "azidomethane"),
    ("CCN=[N+]=[N-]",  "azidoethane"),
    ("CCCN=[N+]=[N-]", "azidopropane"),

    # ── regressions ──────────────────────────────────────────────────────
    ("SC#N",           "thiocyanic acid"),
    ("OC#N",           "cyanic acid"),
    ("CC#N",           "acetonitrile"),
    # acyl azide stays as is
    ("CC(=O)N=[N+]=[N-]", "acetyl azide"),
])
def test_phase293_substitutive_pin(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
