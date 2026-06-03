"""Phase 261: indane-1,3-dione and other fused ring retained names (IUPAC 2013).

  O=C1CC(=O)c2ccccc21 → indane-1,3-dione
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # indane-1,3-dione (diketone fused bicyclic)
    ("O=C1CC(=O)c2ccccc21", "indane-1,3-dione"),
    # regression: related indane compounds
    ("C1CCc2ccccc21",       "indane"),
    ("O=C1CCc2ccccc21",     "indan-1-one"),
    ("O=C1Cc2ccccc2N1",     "indolin-2-one"),
    # regression: phthalimide (N-containing, different from indane-1,3-dione)
    ("O=C1NC(=O)c2ccccc21", "phthalimide"),
])
def test_phase261_indane_dione(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
