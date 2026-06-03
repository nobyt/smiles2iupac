"""Phase 378: Disiloxane and disilazane naming (IUPAC 2013 P-68.4).

Si-O-Si (disiloxane) and Si-NH-Si (disilazane) bonds were previously
undetected, resulting in each Si being named as an independent silane.

Fix: detect Si-O-Si and Si-N-Si pairs in the Si scan loop and emit
disiloxane_org / disilazane_org functional groups covering both Si atoms.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Disiloxane: Si-O-Si
    ("C[Si](C)(O[Si](C)(C)C)C",     "hexamethyldisiloxane"),
    # Disilazane: Si-N(H)-Si
    ("C[Si](C)(N[Si](C)(C)C)C",     "hexamethyldisilazane"),
    # Mixed-alkyl disiloxane
    ("CC[Si](C)(O[Si](C)(C)C)C",    "ethylpentamethyldisiloxane"),
    # Regressions: individual silane/silanol/silyl ether unchanged
    ("[Si](C)(C)(C)C",              "tetramethylsilane"),
    ("[Si](C)(C)(C)O",              "trimethylsilanol"),
    ("[Si](C)(C)(C)OC",             "methoxytrimethylsilane"),
])
def test_phase378_disiloxane_disilazane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
