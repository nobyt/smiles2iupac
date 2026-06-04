"""Phase 401: Uracil and thymine → systematic IUPAC 2013 names.

IUPAC 2013 P-31.1.6.2: retained trivial names uracil and thymine are not
preferred IUPAC names. The preferred names are:
  uracil   → pyrimidine-2,4(1H,3H)-dione
  thymine  → 5-methylpyrimidine-2,4(1H,3H)-dione

Changes:
1. Removed uracil/thymine from the SMILES lookup in __init__.py.
2. Added (True, ("NH","C","C","C","NH","C")) → ("pyrimidine", True) to
   _RETAINED_NAMES so the aromatic tautomer ring is recognised.
3. Modified the first dione block: when is_nh=True, format as
   {base}-{co_locs}({nh_loc1}H,{nh_loc2}H,...)-{dione/trione} using the
   canonical rotation (which already places N at lowest locants).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # uracil (Kekule and aromatic SMILES)
    ("O=C1NC=CC(=O)N1",           "pyrimidine-2,4(1H,3H)-dione"),
    ("O=c1cc[nH]c(=O)[nH]1",      "pyrimidine-2,4(1H,3H)-dione"),
    # thymine (5-methyluracil)
    ("O=C1NC(=O)C(C)=CN1",        "5-methylpyrimidine-2,4(1H,3H)-dione"),
    ("Cc1c[nH]c(=O)[nH]c1=O",     "5-methylpyrimidine-2,4(1H,3H)-dione"),
    # regression: plain pyrimidine unchanged
    ("c1ccncn1",                   "pyrimidine"),
    # regression: mono-ketone pyrimidinones unchanged
    ("O=C1NC=CC=N1",              "pyrimidin-2(1H)-one"),
    ("O=C1NC=NC=C1",              "pyrimidin-4(3H)-one"),
    # regression: cytosine unchanged
    ("NC1=NC(=O)NC=C1",           "4-aminopyrimidin-2(1H)-one"),
    # regression: saturated lactam unchanged
    ("O=C1CCCCN1",                 "piperidin-2-one"),
])
def test_phase401_uracil_thymine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
