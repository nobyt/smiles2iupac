"""Phase 379: Silanediol/silanetriol, tetraalkoxysilane, and disilane.

Three Si-compound fixes:
1. Silanol namer emitted "silanol" regardless of OH count; two OH → "silanediol",
   three OH → "silanetriol".
2. Tetraalkoxysilane (Si with only O-C groups, no direct C-Si) was falling through
   as an unrecognised group.
3. Si-Si bond (disilane) was treated as two independent silanes.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Silanediol / silanetriol
    ("[Si](C)(C)(O)(O)",             "dimethylsilanediol"),
    ("[Si](C)(O)(O)(O)",             "methylsilanetriol"),
    # Tetraalkoxysilane (no direct C-Si)
    ("[Si](OCC)(OCC)(OCC)OCC",       "tetraethoxysilane"),
    ("[Si](OC)(OC)(OC)OC",           "tetramethoxysilane"),
    # Disilane: Si-Si bond (symmetric, fully substituted → no locants)
    ("[Si](C)(C)(C)[Si](C)(C)C",     "hexamethyldisilane"),
    # Regressions
    ("[Si](C)(C)(C)O",               "trimethylsilanol"),
    ("[Si](C)(C)(C)C",               "tetramethylsilane"),
    ("C[Si](C)(O[Si](C)(C)C)C",      "hexamethyldisiloxane"),
    ("[Si](C)(C)(C)OC",              "methoxytrimethylsilane"),
])
def test_phase379_si_compounds(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
