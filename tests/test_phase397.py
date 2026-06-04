"""Phase 397: Pyrimidinone naming (IUPAC 2013 P-31.1.7).

Two fixes:
1. Added ("NH","C","C","C","N","C") → ("pyrimidine", True) to _RETAINED_NAMES.
2. For multi-N indicated-H rings, recompute locants by IUPAC rule:
   minimise (heteroatom_locant_set, co_locant, nh_locant).

Pyrimidin-2(1H)-one: C2=O flanked by N1H and N3 (two N neighbours).
Pyrimidin-4(3H)-one: C4=O flanked by N3H and C5 (one N, one C; N at {1,3} wins).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrimidin-2(1H)-one: C=O flanked by both N atoms
    ("O=C1NC=CC=N1",        "pyrimidin-2(1H)-one"),
    # pyrimidin-4(3H)-one: C=O flanked by NH and C (NH at pos 3)
    ("O=C1NC=NC=C1",        "pyrimidin-4(3H)-one"),
    # cytosine: 4-aminopyrimidin-2(1H)-one
    ("NC1=NC(=O)NC=C1",     "4-aminopyrimidin-2(1H)-one"),
    # regression: plain pyrimidine unchanged
    ("c1ccncn1",             "pyrimidine"),
    # regression: plain pyridine unchanged
    ("c1ccccn1",             "pyridine"),
    # regression: saturated lactam (no is_nh) unchanged
    ("O=C1CCCCN1",           "piperidin-2-one"),
    # regression: 1H-imidazole (no exo C=O) unchanged
    ("c1cnc[nH]1",           "1H-imidazole"),
    # regression: pyridin-2(1H)-one (Phase 396, single-N ring) unchanged
    ("O=C1NC=CC=C1",         "pyridin-2(1H)-one"),
])
def test_phase397_pyrimidinone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
