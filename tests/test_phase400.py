"""Phase 400: 5-membered oxazolone and thiazolone naming (IUPAC 2013 P-31.1.7).

Two changes:
1. Added _RETAINED_NAMES entries for the NH-first tautomers:
     ("NH","C","C","O","C") → ("1,3-oxazole", True)
     ("NH","C","C","S","C") → ("1,3-thiazole", True)
2. Unified the lactam/thiolactam is_nh=True rotation algorithm to a
   priority-aware scheme (O < S < N), replacing the old multi-N vs single-N
   split.  Key: (_o_locs, _s_locs, _n_locs, co_loc, nh_loc) — O gets the
   lowest locant, matching IUPAC P-14.5.2 heteroatom priority rule.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,3-oxazol-2(3H)-one
    ("O=C1NC=CO1",    "1,3-oxazol-2(3H)-one"),
    # 1,3-thiazol-2(3H)-one
    ("O=C1NC=CS1",    "1,3-thiazol-2(3H)-one"),
    # regression: plain parent rings unchanged
    ("c1ncco1",        "oxazole"),
    ("c1cscn1",        "thiazole"),
    # regression: 6-membered indicated-H lactams unchanged
    ("O=C1NC=CC=C1",  "pyridin-2(1H)-one"),
    ("O=C1NC=CC=N1",  "pyrimidin-2(1H)-one"),
    ("O=C1NC=NC=C1",  "pyrimidin-4(3H)-one"),
    ("O=C1NC=CN=C1",  "pyrazin-2(1H)-one"),
    ("O=C1C=CC=NN1",  "pyridazin-3(2H)-one"),
    # regression: 6-membered indicated-H thiolactams unchanged
    ("S=C1NC=CC=C1",  "pyridin-2(1H)-thione"),
    ("S=C1NC=CC=N1",  "pyrimidin-2(1H)-thione"),
    # regression: 1H-imidazole unchanged
    ("c1cnc[nH]1",     "1H-imidazole"),
])
def test_phase400_oxazolone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
