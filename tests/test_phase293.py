"""Phase 293: cyanate/thiocyanate/azide PIN (IUPAC 2013 P-65.3.1).

cyanate/thiocyanate use functional-class names (same as isocyanate):
  COC#N        → methyl cyanate
  CSC#N        → methyl thiocyanate
azide uses substitutive prefix (PIN per P-68.3.1):
  CN=[N+]=[N-] → azidomethane
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── cyanate PIN ───────────────────────────────────────────────────────
    ("COC#N",          "methyl cyanate"),
    ("CCOC#N",         "ethyl cyanate"),
    ("CCCOC#N",        "propyl cyanate"),

    # ── thiocyanate PIN ──────────────────────────────────────────────────
    ("CSC#N",          "methyl thiocyanate"),
    ("CCSC#N",         "ethyl thiocyanate"),
    ("CCCSC#N",        "propyl thiocyanate"),

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
