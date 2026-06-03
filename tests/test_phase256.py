"""Phase 256: carbonyl sulfide and ammonium retained names (IUPAC 2013).

  O=C=S   → carbonyl sulfide  (IUPAC 2013 P-13.1, retained name)
  [NH4+]  → ammonium          (IUPAC 2013 P-73.1.1, retained name)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # carbonyl sulfide
    ("O=C=S",           "carbonyl sulfide"),
    # ammonium
    ("[NH4+]",          "ammonium"),
    # regression: related retained names unchanged
    ("O=C=O",           "carbon dioxide"),
    ("S=C=S",           "carbon disulfide"),
    ("[C-]#[O+]",       "carbon monoxide"),
    ("N",               "ammonia"),
])
def test_phase256_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
