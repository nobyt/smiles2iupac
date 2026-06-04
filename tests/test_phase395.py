"""Phase 395: Thioamide with ring N – {N_locant}-(thioyl){ring_base} (IUPAC 2013 P-65.1.2.2).

Two bugs prevented correct naming of CC(=S)N1CCCCC1:
1. _is_thioamide excluded ring N entirely (and with it, piperidyl thioamides
   were misdetected as thioketone + amine).
2. _name_thioamide returned None when a ring neighbor was detected, even when
   that neighbour was the amide N (not a ring-attached C).

Fix: allow ring N in _is_thioamide when the carbonyl C is not itself in the
ring; restrict the early-return delegation to C-in-ring neighbours only.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thioamide with ring N
    ("CC(=S)N1CCCCC1",      "1-(ethanethioyl)piperidine"),
    ("CC(=S)N1CCCC1",       "1-(ethanethioyl)pyrrolidine"),
    ("CC(=S)N1CCOCC1",      "4-(ethanethioyl)morpholine"),
    # regression: non-ring thioamides unchanged
    ("CC(=S)N",             "ethanethioamide"),
    ("CC(=S)NC",            "N-methylethanethioamide"),
    ("CC(=S)N(C)C",         "N,N-dimethylethanethioamide"),
    # regression: thiolactam (ring C=S) still works (Phase 391)
    ("S=C1CCCCN1",          "piperidine-2-thione"),
    ("S=C1CCCN1",           "pyrrolidine-2-thione"),
    # regression: ring-C neighbour → cycloalkanecarbothioamide still delegates
    ("NC(=S)C1CCCC1",       "cyclopentanecarbothioamide"),
])
def test_phase395_thioamide_ring_n(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
