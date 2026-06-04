"""Phase 394: Cyclic N in sulfonamide → {N_locant}-(sulfonyl){ring} (IUPAC 2013 P-65.3.1).

When the N of a sulfonamide is part of a ring, the code was treating both
ring-C neighbours of N as separate substituents, giving wrong names like
'N-piperidin-2-yl-N-piperidin-6-ylmethanesulfonamide'.

Fix: detect ring N and name as {N-locant}-({sulfonyl_stem}sulfonyl){ring_base}.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ring N in sulfonamide
    ("O=S(=O)(N1CCCCC1)C",    "1-(methanesulfonyl)piperidine"),
    ("O=S(=O)(N1CCCC1)C",     "1-(methanesulfonyl)pyrrolidine"),
    ("O=S(=O)(N1CCOCC1)C",    "4-(methanesulfonyl)morpholine"),
    ("O=S(=O)(N1CCCCC1)CC",   "1-(ethanesulfonyl)piperidine"),
    # regression: non-ring sulfonamides unchanged
    ("CS(=O)(=O)N",            "methanesulfonamide"),
    ("CS(=O)(=O)NC",           "N-methylmethanesulfonamide"),
    ("CS(=O)(=O)N(C)C",        "N,N-dimethylmethanesulfonamide"),
])
def test_phase394_sulfonamide_ring_n(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
