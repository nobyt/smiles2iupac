"""Phase 522: Fix 1H/2H indicated hydrogen assignment for 1,2,3-triazole.

1H-1,2,3-triazole: NH at N1, which is adjacent to one C and one N.
2H-1,2,3-triazole: NH at N2, which is flanked by N1 and N3 (both N).

The lookup table had the two tautomers swapped.
IUPAC 2013 P-31.1.3.4: indicated H locant = position of the NH.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2H form: NH between two N atoms (position 2 in 1,2,3-triazole)
    ("c1cn[nH]n1",  "2H-1,2,3-triazole"),
    # 1H form: NH adjacent to a C atom (position 1 in 1,2,3-triazole)
    ("c1c[nH]nn1",  "1H-1,2,3-triazole"),
    # regression: pyrazole unchanged (NH adjacent to C → 1H)
    ("c1cn[nH]c1",  "1H-pyrazole"),
    ("c1cc[nH]n1",  "1H-pyrazole"),
    # regression: 1H-imidazole unchanged
    ("c1cnc[nH]1",  "1H-imidazole"),
    # regression: 1H-tetrazole unchanged
    ("c1nn[nH]n1",  "1H-tetrazole"),
    # regression: 1H-1,2,4-triazole unchanged
    ("c1nc[nH]n1",  "1H-1,2,4-triazole"),
])
def test_phase522_triazole_indicated_h(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
