"""Phase 396: Pyridinone naming – pyridin-2(1H)-one (IUPAC 2013 P-31.1.7).

For 2-pyridinone (O=C1NC=CC=C1), RDKit treats the ring as aromatic with an
NH, giving canonical sig ('NH','C','C','C','C','C').  Previously this signature
was absent from _RETAINED_NAMES, so name_heterocycle returned None and the
compound fell through to a generic 'oxobenzene' fallback.

Fix:
1. Added ("NH","C","C","C","C","C") → ("pyridine", True) to _RETAINED_NAMES.
2. In the lactam-detection block, when is_nh is True use the indicated-hydrogen
   format: pyridin-{loc}({nh_loc}H)-one instead of 1H-pyridin-{loc}-one.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2-pyridinone (kekulé and aromatic SMILES)
    ("O=C1NC=CC=C1",       "pyridin-2(1H)-one"),
    ("O=c1cccc[nH]1",      "pyridin-2(1H)-one"),
    # regression: plain pyridine unchanged
    ("c1ccccn1",            "pyridine"),
    # regression: saturated lactams use old format (no indicated H needed)
    ("O=C1CCCCN1",          "piperidin-2-one"),
    ("O=C1CCCN1",           "pyrrolidin-2-one"),
    # regression: 1H-imidazole (5-membered NH ring, no exo-ketone) unchanged
    ("c1cnc[nH]1",          "1H-imidazole"),
])
def test_phase396_pyridinone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
