"""Phase 290: methanediimine as PIN for HN=C=NH (IUPAC 2013 P-62.3.1.2).

"carbodiimide" is a retained acceptable name but not the PIN.
The PIN for HN=C=NH is methanediimine (functional replacement nomenclature).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── PIN for the parent compound ───────────────────────────────────────
    ("N=C=N",   "methanediimine"),

    # ── regressions ───────────────────────────────────────────────────────
    ("NNC(=O)N",   "semicarbazide"),
    ("NNC(=S)N",   "thiosemicarbazide"),
    ("NC(=N)N",    "guanidine"),
    ("NC(=O)N",    "urea"),
])
def test_phase290_methanediimine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
