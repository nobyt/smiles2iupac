"""Phase 377: Silyl ester naming and acyloxy contraction fix.

_make_oxy_name incorrectly contracted acyl group names:
  'propanoyl' -> base='propano' + 'oxy' = 'propanooxy' (WRONG)
The fix detects acyl groups (base ends in 'o') and keeps the 'yl':
  'propanoyl' -> 'propanoyl' + 'oxy' = 'propanoyloxy' (CORRECT)

This affected Si-O-C(=O)-R compounds (silyl esters) and any other context
where _make_oxy_name was applied to an acyl group name.

Regression: simple alkyl groups (methyl→methoxy, ethyl→ethoxy) still work.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Silyl esters: Si-O-C(=O)-R
    ("[Si](C)(C)(C)OC(=O)CC",    "trimethylpropanoyloxysilane"),
    ("[Si](C)(C)(C)OC(=O)CCC",   "butanoyloxytrimethylsilane"),
    # Acetoxy still works (acetyl → acetoxy via existing contraction)
    ("[Si](C)(C)(C)OC(=O)C",     "acetoxytrimethylsilane"),
    # Simple silyl ethers unchanged
    ("[Si](C)(C)(C)OC",          "methoxytrimethylsilane"),
    ("[Si](C)(C)(C)OCC",         "ethoxytrimethylsilane"),
    # Regressions: trimethylsilanol unchanged
    ("[Si](C)(C)(C)O",           "trimethylsilanol"),
    # Regressions: trimethylsilane unchanged
    ("[Si](C)(C)(C)C",           "tetramethylsilane"),
])
def test_phase377_silyl_ester_acyloxy(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
