"""Phase 222: selenol and tellurol naming (IUPAC 2013 P-63.4.3).

R-SeH → alkan-N-selenol; R-TeH → alkan-N-tellurol.
No locant for C1 on chains ≤ C2 (analogous to thiol).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # selenol: C1 on C1/C2 chain → no locant
    ("C[SeH]",   "methaneselenol"),
    ("CC[SeH]",  "ethaneselenol"),
    # selenol: longer chain → locant
    ("CCC[SeH]", "propane-1-selenol"),
    # tellurol: analogous
    ("C[TeH]",   "methanetellurol"),
    ("CC[TeH]",  "ethanetellurol"),
    # regression: thiol still works
    ("CS",       "methanethiol"),
    ("CCS",      "ethanethiol"),
])
def test_phase222_selenol_tellurol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
