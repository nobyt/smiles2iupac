"""Phase 380: cyclic sulfate / sulfite ring detection (1,3,2-dioxathiolane / 1,3,2-dioxathiane).

IUPAC 2013 P-31.1.3 retained names for O,O,S heterocycles.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # unsubstituted 5-membered O,O,S ring
    ("C1COSO1",          "1,3,2-dioxathiolane"),
    # 5-membered cyclic sulfite (1 exocyclic =O on S)
    ("O=S1OCCO1",        "2-oxo-1,3,2-dioxathiolane"),
    # 5-membered cyclic sulfate (2 exocyclic =O on S)
    ("O=S1(=O)OCCO1",    "2,2-dioxo-1,3,2-dioxathiolane"),
    # unsubstituted 6-membered O,O,S ring
    ("C1CCOSO1",         "1,3,2-dioxathiane"),
    # 6-membered cyclic sulfite
    ("O=S1OCCCO1",       "2-oxo-1,3,2-dioxathiane"),
    # 6-membered cyclic sulfate
    ("O=S1(=O)OCCCO1",   "2,2-dioxo-1,3,2-dioxathiane"),
])
def test_phase380_cyclic_sulfate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
