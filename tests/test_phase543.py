"""Phase 543: S-aryl thioate missing parentheses around locant-bearing substituent.

_name_thioester returned "S-pyridin-2-yl ethanethioate" — the substituent name
"pyridin-2-yl" contains a locant and requires enclosing parentheses per IUPAC.
Simple names (phenyl, methyl, ethyl) with no locants are unchanged.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # S-heteroaryl thioates: need parens around locant-bearing name
    ("O=C(Sc1ccccn1)C",   "S-(pyridin-2-yl) ethanethioate"),
    ("O=C(Sc1cccnc1)C",   "S-(pyridin-3-yl) ethanethioate"),
    ("O=C(Sc1ccncc1)C",   "S-(pyridin-4-yl) ethanethioate"),
    ("O=C(Sc1cccs1)C",    "S-(thiophen-2-yl) ethanethioate"),
    ("O=C(Sc1ccccn1)CC",  "S-(pyridin-2-yl) propanethioate"),
    # Regression: simple substituents unchanged (no digits → no parens)
    ("O=C(Sc1ccccc1)C",   "S-phenyl ethanethioate"),
    ("O=C(SC)CC",          "S-methyl propanethioate"),
    ("O=C(SCC)C",          "S-ethyl ethanethioate"),
])
def test_phase543_s_aryl_thioate_parens(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
